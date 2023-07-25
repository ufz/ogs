/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on September 22, 2022, 12:29 PM
 */

#include <tclap/CmdLine.h>

#include <algorithm>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "BaseLib/FileTools.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Location.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "MeshToolsLib/IntegrationPointDataTools.h"

struct MeshEntityMapInfo
{
    std::size_t const partition_id;
    std::size_t const original_id;
};

template <typename T>
bool createPropertyVector(
    MeshLib::Mesh& merged_mesh,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& partitioned_meshes,
    MeshLib::PropertyVector<T> const* const pv,
    MeshLib::Properties const& properties,
    std::vector<MeshEntityMapInfo> const& merged_node_map,
    std::vector<MeshEntityMapInfo> const& merged_element_map)
{
    if (pv == nullptr)
    {
        return false;
    }

    if (pv->getPropertyName() == "vtkGhostType")
    {
        // Do nothing
        return true;
    }

    auto const item_type = pv->getMeshItemType();

    auto const pv_name = pv->getPropertyName();

    auto const pv_num_components = pv->getNumberOfGlobalComponents();

    if (pv_name == "OGS_VERSION" || pv_name == "IntegrationPointMetaData")
    {
        auto new_pv = MeshLib::getOrCreateMeshProperty<T>(
            merged_mesh, pv_name, item_type, pv_num_components);
        new_pv->resize(pv->size());

        std::copy(pv->begin(), pv->end(), new_pv->begin());

        return true;
    }

    std::vector<MeshLib::PropertyVector<T>*> partition_property_vectors;
    partition_property_vectors.reserve(partitioned_meshes.size());
    for (auto const& mesh : partitioned_meshes)
    {
        partition_property_vectors.emplace_back(
            mesh->getProperties().getPropertyVector<T>(pv_name, item_type,
                                                       pv_num_components));
    }

    auto createNewCellOrNodePropertyVector =
        [&](std::vector<MeshEntityMapInfo> const& mesh_entity_map)
    {
        auto new_pv = MeshLib::getOrCreateMeshProperty<T>(
            merged_mesh, pv_name, item_type, pv_num_components);
        std::size_t counter = 0;
        for (auto const& entity_info : mesh_entity_map)
        {
            auto const& partition_pv =
                partition_property_vectors[entity_info.partition_id];
            for (int i_com = 0; i_com < pv_num_components; i_com++)
            {
                (*new_pv)[counter * pv_num_components + i_com] =
                    (*partition_pv)[entity_info.original_id *
                                        pv_num_components +
                                    i_com];
            }
            counter++;
        }
    };

    if (item_type == MeshLib::MeshItemType::Node)
    {
        createNewCellOrNodePropertyVector(merged_node_map);
        return true;
    }

    if (item_type == MeshLib::MeshItemType::Cell)
    {
        createNewCellOrNodePropertyVector(merged_element_map);
        return true;
    }

    if (item_type == MeshLib::MeshItemType::IntegrationPoint)
    {
        std::vector<std::vector<std::size_t>> partition_element_offsets;
        partition_element_offsets.reserve(partitioned_meshes.size());
        for (auto const& mesh : partitioned_meshes)
        {
            partition_element_offsets.emplace_back(
                MeshToolsLib::getIntegrationPointDataOffsetsOfMeshElements(
                    mesh->getElements(), *pv, properties));
        }

        auto new_pv = MeshLib::getOrCreateMeshProperty<T>(
            merged_mesh, pv_name, item_type, pv_num_components);

        // Count the integration points
        std::size_t counter = 0;
        auto const ip_meta_data =
            MeshLib::getIntegrationPointMetaData(properties, pv_name);

        for (auto const element : merged_mesh.getElements())
        {
            int const number_of_integration_points =
                MeshToolsLib::getNumberOfElementIntegrationPoints(ip_meta_data,
                                                                  *element);
            counter += number_of_integration_points;
        }
        new_pv->resize(counter * pv_num_components);

        auto const global_ip_offsets =
            MeshToolsLib::getIntegrationPointDataOffsetsOfMeshElements(
                merged_mesh.getElements(), *pv, properties);

        std::size_t element_counter = 0;
        for (auto const& element_info : merged_element_map)
        {
            MeshLib::PropertyVector<T> const& partition_pv =
                *(partition_property_vectors[element_info.partition_id]);

            auto const& offsets =
                partition_element_offsets[element_info.partition_id];

            int const begin_pos = offsets[element_info.original_id];
            int const end_pos = offsets[element_info.original_id + 1];

            std::copy(partition_pv.begin() + begin_pos,
                      partition_pv.begin() + end_pos,
                      new_pv->begin() + global_ip_offsets[element_counter]);

            element_counter++;
        }
    }

    return true;
}

std::vector<std::string> readVtuFileNames(std::string const& pvtu_file_name)
{
    std::ifstream ins(pvtu_file_name);

    if (!ins)
    {
        OGS_FATAL("Could not open pvtu file {:s}.", pvtu_file_name);
    }

    using boost::property_tree::ptree;
    ptree pt;
    read_xml(ins, pt);

    auto root = pt.get_child("VTKFile");

    std::vector<std::string> vtu_file_names;

    std::string file_path = BaseLib::extractPath(pvtu_file_name);

    for (ptree::value_type const& v : root.get_child("PUnstructuredGrid"))
    {
        if (v.first == "Piece")
        {
            vtu_file_names.push_back(BaseLib::joinPaths(
                file_path,
                // only gets the vtu file name:
                std::filesystem::path(v.second.get("<xmlattr>.Source", ""))
                    .filename()
                    .string()));
        }
    }

    if (vtu_file_names.empty())
    {
        OGS_FATAL("PVTU file {:s} does not contain any vtu piece",
                  pvtu_file_name);
    }

    return vtu_file_names;
}

std::tuple<std::vector<MeshLib::Element*>, std::vector<MeshEntityMapInfo>>
getRegularElements(std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    std::vector<MeshLib::Element*> regular_elements;

    std::size_t partition_counter = 0;
    std::vector<MeshEntityMapInfo> merged_element_map;
    for (auto const& mesh : meshes)
    {
        MeshLib::Properties const& properties = mesh->getProperties();

        std::vector<unsigned char> const& ghost_id_vector =
            *(properties.getPropertyVector<unsigned char>("vtkGhostType"));

        auto const& mesh_elements = mesh->getElements();

        auto const last_element_id_of_previous_partition =
            regular_elements.size();

        std::copy_if(mesh_elements.begin(), mesh_elements.end(),
                     std::back_inserter(regular_elements),
                     [&ghost_id_vector](auto const element)
                     { return ghost_id_vector[element->getID()] == 0; });

        for (auto element_id = last_element_id_of_previous_partition;
             element_id < regular_elements.size();
             element_id++)
        {
            merged_element_map.push_back(
                {partition_counter, regular_elements[element_id]->getID()});
        }

        partition_counter++;
    }

    return {regular_elements, merged_element_map};
}

std::tuple<std::vector<MeshLib::Node*>, std::vector<MeshEntityMapInfo>>
getNodesOfRegularElements(
    std::vector<MeshLib::Element*> const& regular_elements,
    std::vector<MeshEntityMapInfo> const& merged_element_map)
{
    std::vector<MeshLib::Node*> merged_nodes;
    std::vector<MeshEntityMapInfo> merged_node_map;
    std::unordered_set<const MeshLib::Node*> node_status;

    std::size_t element_counter = 0;
    for (auto& e : regular_elements)
    {
        for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
        {
            auto* const node = e->getNode(i);

            if (node_status.contains(node))

            {
                continue;
            }

            // TODO To use std::execution::par if Mac compiler supports
            // parallel algorithms.
            auto const found_node =
                std::find_if(std::begin(merged_nodes),
                             std::end(merged_nodes),
                             [&node](auto const merged_node) {
                                 return (node->asEigenVector3d() ==
                                         merged_node->asEigenVector3d());
                             });

            if (found_node != std::end(merged_nodes))

            {
                e->setNode(i, *found_node);
                continue;
            }

            node_status.insert(node);
            merged_nodes.push_back(node);
            merged_node_map.push_back(
                {merged_element_map[element_counter].partition_id,
                 node->getID()});
        }
        element_counter++;
    }

    return {merged_nodes, merged_node_map};
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "This tool merges VTU files of PVTU into one single VTU file. Apart "
        "from the mesh data, all property data are merged as well"
        "\n\nOpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2023, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "the output mesh (*.vtu)", true, "", "output.vtu");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "the partitioned input mesh (*.pvtu)", true, "",
        "input.pvtu");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    if (BaseLib::getFileExtension(input_arg.getValue()) != ".pvtu")
    {
        OGS_FATAL("The extension of input file name {:s} is not \"pvtu\"",
                  input_arg.getValue());
    }
    if (BaseLib::getFileExtension(output_arg.getValue()) != ".vtu")
    {
        OGS_FATAL("The extension of output file name {:s} is not \"vtu\"",
                  output_arg.getValue());
    }

    auto const vtu_file_names = readVtuFileNames(input_arg.getValue());

    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    meshes.reserve(vtu_file_names.size());

    for (auto const& file_name : vtu_file_names)
    {
        auto mesh = std::unique_ptr<MeshLib::Mesh>(
            MeshLib::IO::VtuInterface::readVTUFile(file_name));

        MeshLib::Properties const& properties = mesh->getProperties();

        if (!properties.existsPropertyVector<unsigned char>("vtkGhostType"))
        {
            OGS_FATAL(
                "Property vector vtkGhostType does not exist in mesh {:s}.",
                file_name);
        }

        meshes.emplace_back(std::move(mesh));
    }

    // If structured binding is used for the returned tuple, Mac compiler gives
    // an error in reference to local binding in calling applyToPropertyVectors.
    std::vector<MeshEntityMapInfo> merged_element_map;
    std::vector<MeshLib::Element*> regular_elements;
    std::tie(regular_elements, merged_element_map) = getRegularElements(meshes);

    std::vector<MeshEntityMapInfo> merged_node_map;
    std::vector<MeshLib::Node*> merged_nodes;
    std::tie(merged_nodes, merged_node_map) =
        getNodesOfRegularElements(regular_elements, merged_element_map);

    // The Node pointers of 'merged_nodes' and Element pointers of
    // 'regular_elements' are shared with 'meshes', the partitioned meshes.
    MeshLib::Mesh merged_mesh =
        MeshLib::Mesh("pvtu_merged_mesh", merged_nodes, regular_elements);

    auto const& properties = meshes[0]->getProperties();

    applyToPropertyVectors(
        properties,
        [&](auto type, auto const& property)
        {
            return createPropertyVector<decltype(type)>(
                merged_mesh, meshes,
                dynamic_cast<MeshLib::PropertyVector<decltype(type)> const*>(
                    property),
                properties, merged_node_map, merged_element_map);
        });

    MeshLib::IO::VtuInterface writer(&merged_mesh);

    auto const result = writer.writeToFile(output_arg.getValue());
    if (!result)
    {
        ERR("Could not write mesh to '{:s}'.", output_arg.getValue());
        return EXIT_FAILURE;
    }

    // Since the Node pointers of 'merged_nodes' and Element pointers of
    // 'regular_elements' are held by 'meshes', the partitioned meshes, the
    // memory by these pointers are released by 'meshes' automatically.
    // Therefore, only node vector and element vector of merged_mesh should be
    // cleaned.
    merged_mesh.shallowClean();

    return EXIT_SUCCESS;
}
