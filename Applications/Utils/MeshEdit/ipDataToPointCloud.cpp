/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <unordered_map>

#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixCache.h"

template <typename MeshElement>
std::vector<std::array<double, 3>> interpolateElementCoords(
    MeshLib::Element const& e, NumLib::ShapeMatrixCache const& sm_cache)
{
    constexpr int dim = MeshElement::dimension;
    using ShapeFunction =
        NumLib::ShapeMatrixCache::ShapeFunctionHigherOrder<MeshElement>;
    using ShapeMatrices = EigenFixedShapeMatrixPolicy<ShapeFunction, dim>;

    auto const& Ns = sm_cache.NsHigherOrder<MeshElement>();
    std::vector<std::array<double, 3>> coords;
    coords.reserve(Ns.size());

    for (auto const& N : Ns)
    {
        coords.push_back(
            NumLib::interpolateCoordinates<ShapeFunction, ShapeMatrices>(e, N));
    }

    return coords;
}

// A template metafunction "returning" the cell type of a mesh element
template <typename Element>
struct CellTypeOfTemplateElement;

template <typename ElementRule>
struct CellTypeOfTemplateElement<MeshLib::TemplateElement<ElementRule>>
{
    static constexpr MeshLib::CellType value = ElementRule::cell_type;
};

std::unordered_map<MeshLib::CellType, std::vector<std::array<double, 3>> (*)(
                                          MeshLib::Element const&,
                                          NumLib::ShapeMatrixCache const&)>
createElementCoordInterpolatorsForAllElementTypes()
{
    std::unordered_map<MeshLib::CellType,
                       std::vector<std::array<double, 3>> (*)(
                           MeshLib::Element const&,
                           NumLib::ShapeMatrixCache const&)>
        map_cell_type_to_element_coords_interpolator;

    map_cell_type_to_element_coords_interpolator.reserve(
        boost::mp11::mp_size<NumLib::AllElementTraitsLagrange>::value);

    boost::mp11::mp_for_each<NumLib::AllElementTraitsLagrange>(
        [&map_cell_type_to_element_coords_interpolator]<typename ETL>(ETL)
        {
            using MeshElement = typename ETL::Element;
            auto constexpr cell_type =
                CellTypeOfTemplateElement<MeshElement>::value;

            auto const [it, newly_inserted] =
                map_cell_type_to_element_coords_interpolator.emplace(
                    cell_type, interpolateElementCoords<MeshElement>);

            if (!newly_inserted)
            {
                OGS_FATAL(
                    "Coordinate interpolator for cell type {} is already "
                    "present",
                    MeshLib::CellType2String(cell_type));
            }
        });

    return map_cell_type_to_element_coords_interpolator;
}

std::vector<MeshLib::Node*> computePointCloudNodes(
    MeshLib::Mesh const& mesh, unsigned const integration_order)
{
    NumLib::ShapeMatrixCache sm_cache{integration_order};
    auto const map_cell_type_to_element_coords_interpolator =
        createElementCoordInterpolatorsForAllElementTypes();

    std::vector<MeshLib::Node*> nodes;

    for (auto const* element : mesh.getElements())
    {
        auto const cell_type = element->getCellType();

        auto const it =
            map_cell_type_to_element_coords_interpolator.find(cell_type);
        if (it == map_cell_type_to_element_coords_interpolator.end())
        {
            OGS_FATAL("Unsupported cell type {} for element #{}",
                      MeshLib::CellType2String(cell_type), element->getID());
        }

        auto const& element_coords_interpolator = it->second;
        auto const coords = element_coords_interpolator(*element, sm_cache);

        for (auto& cs : coords)
        {
            nodes.push_back(new MeshLib::Node(cs, nodes.size()));
        }
    }

    return nodes;
}

unsigned determineIntegrationOrder(MeshLib::Mesh const& mesh)
{
    auto const& properties = mesh.getProperties();

    std::optional<unsigned> integration_order;

    for (auto const& [name, property] : properties)
    {
        if (property->getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        if (name == "IntegrationPointMetaData" || name == "OGS_VERSION")
        {
            continue;
        }

        auto const order =
            MeshLib::getIntegrationPointMetaData(properties, std::string(name))
                .integration_order;

        if (!integration_order)
        {
            integration_order = order;
        }
        else if (integration_order != order)
        {
            OGS_FATAL("Integration orders differ: {} != {}", *integration_order,
                      order);
        }
    }

    if (!integration_order)
    {
        OGS_FATAL("Integration order could not be determined.");
    }

    return *integration_order;
}

void copyDoubleValuedFieldDataToPointCloud(MeshLib::Properties const& props_in,
                                           MeshLib::Properties& props_out,
                                           std::size_t const num_ips)
{
    for (auto const& [prop_name, prop_in] : props_in)
    {
        auto const mesh_item_type = prop_in->getMeshItemType();

        if (mesh_item_type != MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const num_comp = prop_in->getNumberOfGlobalComponents();
        auto const* prop_in_double =
            dynamic_cast<MeshLib::PropertyVector<double> const*>(prop_in);

        if (prop_in_double == nullptr)
        {
            INFO(
                "Ignoring non-double-valued property '{}' with {} components "
                "on {}",
                prop_name, num_comp, toString(mesh_item_type));
            continue;
        }

        DBUG("Converting property '{}' with {} components on {}", prop_name,
             num_comp, toString(mesh_item_type));

        if (props_out.existsPropertyVector<double>(prop_name))
        {
            OGS_FATAL(
                "A property vector with name '{}' already exists. Not "
                "adding it again",
                prop_name);
        }

        if (auto const num_ips_actual = prop_in_double->getNumberOfTuples();
            num_ips_actual != num_ips)
        {
            WARN(
                "Property vector '{}' has {} tuples, which differs from "
                "the number of integration points ({}). Skipping this "
                "property.",
                prop_name, num_ips_actual, num_ips);
        }

        auto* prop_out = props_out.createNewPropertyVector<double>(
            prop_name, MeshLib::MeshItemType::Node, num_comp);

        static_cast<std::vector<double>&>(*prop_out) = *prop_in_double;
    }
}

int main(int argc, char** argv)
{
    // TODO future additions to this tool might include:
    // -C --copy-cell-data
    // -N --interpolate-node-data
    // -I --add-cell-ids
    // -O --integration-order
    // --natural-coords add natural coordinates of integration points, or better
    //     not?
    // --no-data
    // --linear-shape-functions={pressure, ...} for the interpolation of nodal
    //     data
    // PVD support
    // handle other data than only double data

    TCLAP::CmdLine cmd(
        "Creates a point cloud mesh for the integration point data existing on "
        "a given input mesh.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> arg_out_file(
        "o", "output-file", "the output mesh (point cloud, VTU file)", true, "",
        "path");
    cmd.add(arg_out_file);
    TCLAP::ValueArg<std::string> arg_in_file(
        "i", "input-file", "the input mesh (VTU file)", true, "", "path");
    cmd.add(arg_in_file);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh const> mesh_in(
        MeshLib::IO::readMeshFromFile(arg_in_file.getValue()));

    auto const integration_order = determineIntegrationOrder(*mesh_in);
    auto nodes = computePointCloudNodes(*mesh_in, integration_order);

    MeshLib::Mesh point_cloud("point_cloud", std::move(nodes), {});

    copyDoubleValuedFieldDataToPointCloud(mesh_in->getProperties(),
                                          point_cloud.getProperties(),
                                          point_cloud.getNumberOfNodes());

    MeshLib::IO::writeMeshToFile(point_cloud, arg_out_file.getValue());

    return EXIT_SUCCESS;
}
