/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#include <fstream>

#include "BaseLib/MPI.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

/// Returns a vector of values where each value is associated with a
/// particular node. Since a node is part of elements, it is possible to
/// assign an area per element to this node. Each value of the return vector
/// is the sum of the assigned area (per element) multiplied by the property
/// value of the element.
/// @param mesh a surface mesh containing a property \c prop_name assigned
/// to cells
/// @param prop_name name of the cell based property within the \c mesh
/// @return vector of integration values associated to the surface mesh nodes
std::vector<double> getSurfaceIntegratedValuesForNodes(
    const MeshLib::Mesh& mesh, std::string const& prop_name)
{
    if (mesh.getDimension() != 2)
    {
        ERR("Error in "
            "MeshSurfaceExtraction::getSurfaceIntegratedValuesForNodes() - "
            "Given mesh is no surface mesh (dimension != 2).");
        return std::vector<double>();
    }

    if (!mesh.getProperties().existsPropertyVector<double>(prop_name))
    {
        ERR("Need element property, but the property '{:s}' is not available.",
            prop_name);
        return std::vector<double>();
    }
    auto const* const elem_pv = mesh.getProperties().getPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell, 1);

    std::vector<double> integrated_node_area_vec;
    double total_area(0);

    for (auto const* node : mesh.getNodes())
    {
        double integrated_node_area(0);
        for (auto const& connected_elem :
             mesh.getElementsConnectedToNode(*node))
        {
            double const area = connected_elem->getContent() /
                                connected_elem->getNumberOfBaseNodes();
            integrated_node_area += area * (*elem_pv)[connected_elem->getID()];
            total_area += area;
        }

        integrated_node_area_vec.push_back(integrated_node_area);
    }

    INFO("Total surface area: {:g}", total_area);

    return integrated_node_area_vec;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Integrates the given element property and outputs an OGS-5 direct "
        "Neumann boundary condition. The mesh has to contain a property "
        "'bulk_node_ids' that stores the original subsurface mesh node ids. "
        "Such surface meshes can be created using the OGS-6 tool "
        "ExtractSurface.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> in_mesh(
        "i",
        "in-mesh",
        "the surface mesh that has an element property for the Neumann "
        "boundary condition",
        true,
        "",
        "filename for surface mesh input");
    cmd.add(in_mesh);

    TCLAP::ValueArg<std::string> property_in_arg(
        "p",
        "property-in-name",
        "name of an element property used for the computation of the Neumann "
        "boundary condition",
        true,
        "",
        "string (property name)");
    cmd.add(property_in_arg);

    TCLAP::ValueArg<std::string> property_out_arg(
        "",
        "property-out-name",
        "name of the node based property used for the output of the Neumann "
        "boundary condition",
        true,
        "",
        "string (property name)");
    cmd.add(property_out_arg);

    TCLAP::ValueArg<std::string> result_file(
        "o",
        "result-out",
        "the file name the result will be written to ",
        true,
        "",
        "output file name");
    cmd.add(result_file);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    // read surface mesh
    std::unique_ptr<MeshLib::Mesh> surface_mesh(
        MeshLib::IO::readMeshFromFile(in_mesh.getValue()));

    auto const* const node_id_pv =
        [&]() -> MeshLib::PropertyVector<std::size_t>*
    {
        try
        {
            return surface_mesh->getProperties().getPropertyVector<std::size_t>(
                MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
                MeshLib::MeshItemType::Node, 1);
        }
        catch (std::runtime_error const& e)
        {
            WARN("{:s}", e.what());
            return nullptr;
        }
    }();
    if (!node_id_pv)
    {
        return EXIT_FAILURE;
    }

    std::vector<double> integrated_values = getSurfaceIntegratedValuesForNodes(
        *surface_mesh, property_in_arg.getValue());
    std::vector<std::pair<std::size_t, double>> direct_values;
    direct_values.reserve(surface_mesh->getNumberOfNodes());

    for (auto const* node : surface_mesh->getNodes())
    {
        auto const id(node->getID());
        auto const subsurface_node_id((*node_id_pv)[id]);
        auto const val(integrated_values[id]);
        direct_values.emplace_back(subsurface_node_id, val);
    }

    auto* const pv =
        surface_mesh->getProperties().createNewPropertyVector<double>(
            property_out_arg.getValue(), MeshLib::MeshItemType::Node,
            surface_mesh->getNodes().size(), 1);
    for (std::size_t k(0); k < surface_mesh->getNodes().size(); ++k)
    {
        (*pv)[k] = direct_values[k].second;
    }

    MeshLib::IO::writeMeshToFile(*surface_mesh, result_file.getValue());

    std::ofstream result_out(result_file.getValue() + ".txt");
    result_out.precision(std::numeric_limits<double>::max_digits10);
    for (auto const& p : direct_values)
    {
        result_out << p.first << " " << p.second << "\n";
    }

    return EXIT_SUCCESS;
}
