/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <memory>
#include <string>

#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

std::vector<std::size_t> generateBulkIDsReverseMapping(
    std::vector<std::size_t> const& bulk_ids)
{
    std::vector<std::size_t> parallel_to_sequential_node_ids(bulk_ids.size());

    std::size_t pos = 0;
    for (auto const id : bulk_ids)
    {
        parallel_to_sequential_node_ids[id] = pos;
        pos++;
    }
    return parallel_to_sequential_node_ids;
}

template <typename T>
bool reorderProperty(MeshLib::PropertyVector<T> const& original_pv,
                     std::vector<std::size_t> const& bulk_ids,
                     MeshLib::PropertyVector<T>& reordered_pv)
{
    if (original_pv.getNumberOfGlobalComponents() !=
        reordered_pv.getNumberOfGlobalComponents())
    {
        OGS_FATAL(
            "Could not reorder property {} since the number of components "
            "isn't equal ({} != {}).",
            original_pv.getPropertyName(),
            original_pv.getNumberOfGlobalComponents(),
            reordered_pv.getNumberOfGlobalComponents());
    }

    auto const number_of_components = original_pv.getNumberOfGlobalComponents();
    std::size_t pos = 0;
    for (auto const& id : bulk_ids)
    {
        for (int c = 0; c < number_of_components; c++)
        {
            reordered_pv.getComponent(id, c) = original_pv.getComponent(pos, c);
        }
        pos++;
    }
    return true;
}

template <typename T>
bool reorderProperty(MeshLib::Properties const& original_properties,
                     std::vector<std::size_t> const& bulk_ids,
                     std::string const& property_name,
                     MeshLib::Properties& reordered_properties)
{
    if (!reordered_properties.existsPropertyVector<T>(property_name))
    {
        return false;
    }
    auto const& original_property =
        *original_properties.getPropertyVector<T>(property_name);
    auto& reordered_property =
        *reordered_properties.getPropertyVector<T>(property_name);
    return reorderProperty(original_property, bulk_ids, reordered_property);
}

void reorderProperties(MeshLib::Properties const& original_properties,
                       std::vector<std::size_t> const& bulk_ids,
                       std::vector<std::string> const& property_names,
                       MeshLib::Properties& reordered_properties)
{
    for (auto const& property_name : property_names)
    {
        bool const success =
            reorderProperty<double>(original_properties, bulk_ids,
                                    property_name, reordered_properties) ||
            reorderProperty<float>(original_properties, bulk_ids, property_name,
                                   reordered_properties) ||
            reorderProperty<int>(original_properties, bulk_ids, property_name,
                                 reordered_properties) ||
            reorderProperty<long>(original_properties, bulk_ids, property_name,
                                  reordered_properties) ||
            reorderProperty<unsigned>(original_properties, bulk_ids,
                                      property_name, reordered_properties) ||
            reorderProperty<unsigned long>(original_properties, bulk_ids,
                                           property_name,
                                           reordered_properties) ||
            reorderProperty<std::size_t>(original_properties, bulk_ids,
                                         property_name, reordered_properties) ||
            reorderProperty<char>(original_properties, bulk_ids, property_name,
                                  reordered_properties) ||
            reorderProperty<unsigned char>(original_properties, bulk_ids,
                                           property_name, reordered_properties);
        if (!success)
        {
            OGS_FATAL("Could not reorder PropertyVector '{}'.", property_name);
        }
    }
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reorder mesh nodes and cells.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_in_arg(
        "i", "input", "name of the input mesh file", true, "", "string");
    cmd.add(mesh_in_arg);
    TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output", "name of the output mesh file", true, "", "string");
    cmd.add(mesh_out_arg);

    cmd.parse(argc, argv);

    const std::string filename(mesh_in_arg.getValue());

    // read the mesh file
    auto const mesh =
        std::unique_ptr<MeshLib::Mesh>(MeshLib::IO::readMeshFromFile(filename));
    if (!mesh)
    {
        return EXIT_FAILURE;
    }

    auto const* bulk_node_ids = MeshLib::bulkNodeIDs(*mesh);
    if (!bulk_node_ids)
    {
        OGS_FATAL("Property / data array '{}' has not been found in the mesh.",
                  MeshLib::getBulkIDString(MeshLib::MeshItemType::Node));
    }

    auto const& nodes = mesh->getNodes();
    auto const node_ids_reverse_mapping(
        generateBulkIDsReverseMapping(*bulk_node_ids));

    std::vector<MeshLib::Node*> reordered_nodes(nodes.size());
    std::size_t pos = 0;
    for (auto const id : node_ids_reverse_mapping)
    {
        auto const& n = *(nodes[id]);
        reordered_nodes[pos] = new MeshLib::Node(n[0], n[1], n[2], pos);
        pos++;
    }

    auto const bulk_element_ids_string =
        MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell);
    if (!mesh->getProperties().existsPropertyVector<std::size_t>(
            bulk_element_ids_string))
    {
        OGS_FATAL("Property / data array '{}' has not been found in the mesh.",
                  bulk_element_ids_string);
    }
    auto const& bulk_element_ids =
        *mesh->getProperties().getPropertyVector<std::size_t>(
            bulk_element_ids_string);

    std::vector<std::size_t> element_ids_reverse_mapping(
        generateBulkIDsReverseMapping(bulk_element_ids));

    auto const& elements = mesh->getElements();
    std::vector<MeshLib::Element*> reordered_elements(elements.size());
    pos = 0;
    for (auto const id : element_ids_reverse_mapping)
    {
        auto const& e = *elements[id];
        auto* reordered_element =
            e.clone();  // use clone to have the same element type
        // reset the nodes to the reordered nodes
        auto const number_of_nodes = reordered_element->getNumberOfNodes();
        for (std::size_t node_number = 0; node_number < number_of_nodes;
             node_number++)
        {
            reordered_element->setNode(
                node_number,
                reordered_nodes[(
                    *bulk_node_ids)[e.getNode(node_number)->getID()]]);
        }
        reordered_elements[pos] = reordered_element;
        pos++;
    }
    auto const& original_properties = mesh->getProperties();
    // reordering of PropertyVector with following MeshItemTypes isn't supported
    std::vector<MeshLib::MeshItemType> const exclude_property_vectors{
        MeshLib::MeshItemType::Edge, MeshLib::MeshItemType::Face,
        MeshLib::MeshItemType::IntegrationPoint};
    // Create reordered mesh - only the PropertyVectors for which the reordering
    // is implemented are copied
    MeshLib::Mesh reordered_mesh{
        "reordered_mesh", reordered_nodes, reordered_elements,
        false /* compute_element_neighbors */,
        original_properties.excludeCopyProperties(exclude_property_vectors)};
    auto& properties = reordered_mesh.getProperties();

    // node based properties
    auto const node_property_names =
        mesh->getProperties().getPropertyVectorNames(
            MeshLib::MeshItemType::Node);
    reorderProperties(original_properties, *bulk_node_ids, node_property_names,
                      properties);

    // element based properties
    auto const element_property_names =
        mesh->getProperties().getPropertyVectorNames(
            MeshLib::MeshItemType::Cell);
    reorderProperties(original_properties, bulk_element_ids,
                      element_property_names, properties);
    MeshLib::IO::writeMeshToFile(reordered_mesh, mesh_out_arg.getValue());

    auto const number_of_properties =
        mesh->getProperties().getPropertyVectorNames().size();
    if (node_property_names.size() + element_property_names.size() <
        number_of_properties)
    {
        WARN(
            "Properties that are not assigened to nodes or elements are not "
            "transferred to the reordered mesh.");
    }

    return EXIT_SUCCESS;
}
