/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <map>
#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshNodeSearcher.h"

namespace
{
/// Find one-to-one mapping of subdomain nodes to the bulk mesh nodes and
/// returns the map of ids.
std::vector<std::size_t> identifySubdomainMeshNodes(
    MeshLib::Mesh const& subdomain_mesh,
    MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher)
{
    // Convert nodes pointers needed for the mesh_node_searcher algorithm.
    auto const& nodes = subdomain_mesh.getNodes();
    std::vector<MathLib::Point3dWithID*> subdomain_points{begin(nodes),
                                                          end(nodes)};

    auto const& bulk_node_ids =
        mesh_node_searcher.getMeshNodeIDs(subdomain_points);

    if (bulk_node_ids.size() != subdomain_mesh.getNumberOfNodes())
    {
        OGS_FATAL(
            "Expected to find exactly one node in the bulk mesh for each node "
            "of the subdomain; Found {:d} nodes in the bulk mesh out of {:d} "
            "nodes in the subdomain.",
            bulk_node_ids.size(), subdomain_mesh.getNumberOfNodes());
    }

    return bulk_node_ids;
}

/// Find all elements which consist of all of the given node ids.
/// \note It is recommended to include only base nodes of elements into the
/// search \c node_ids.
std::vector<std::size_t> findElementsInMesh(
    MeshLib::Mesh const& mesh, std::vector<std::size_t> const& node_ids)
{
    //
    // Collect all element ids for all nodes.
    //
    std::vector<std::size_t> common_element_ids;
    // Every node is connected to at least one element.
    auto const nnodes = node_ids.size();
    common_element_ids.reserve(nnodes);

    for (auto const node_id : node_ids)
    {
        auto const& connected_elements = mesh.getNode(node_id)->getElements();
        std::transform(
            begin(connected_elements), end(connected_elements),
            back_inserter(common_element_ids),
            [](MeshLib::Element const* const e) { return e->getID(); });
    }

    //
    // Count how often an element is shared by all nodes.
    //
    std::map<std::size_t, int> element_counts;
    for (auto const element_id : common_element_ids)
    {
        element_counts[element_id]++;
    }

    //
    // Elements which are shared by as many nodes as the input nodes are the
    // desired elements.
    //
    std::vector<std::size_t> element_ids;
    for (auto const& pair : element_counts)
    {
        if (pair.second == static_cast<int>(nnodes))
        {
            element_ids.push_back(pair.first);
        }
    }

    return element_ids;
}

/// Tries to find all elements' ids in the bulk mesh.
/// The subdomain elements' nodes are mapped to the bulk mesh and are then used
/// to identify the bulk mesh elements.
std::vector<std::vector<std::size_t>> identifySubdomainMeshElements(
    MeshLib::Mesh const& subdomain_mesh, MeshLib::Mesh const& bulk_mesh)
{
    auto& properties = subdomain_mesh.getProperties();
    auto const& bulk_node_ids = *properties.getPropertyVector<std::size_t>(
        "bulk_node_ids", MeshLib::MeshItemType::Node, 1);

    // Allocate space for all elements for random insertion.
    std::vector<std::vector<std::size_t>> bulk_element_ids_map(
        subdomain_mesh.getNumberOfElements());

    for (auto* const e : subdomain_mesh.getElements())
    {
        std::vector<std::size_t> element_node_ids(e->getNumberOfBaseNodes());
        for (unsigned n = 0; n < e->getNumberOfBaseNodes(); ++n)
        {
            element_node_ids[n] = e->getNodeIndex(n);
        }
        std::vector<std::size_t> element_node_ids_bulk(
            e->getNumberOfBaseNodes());
        std::transform(begin(element_node_ids), end(element_node_ids),
                       begin(element_node_ids_bulk),
                       [&bulk_node_ids](std::size_t const id) {
                           return bulk_node_ids[id];
                       });
        std::vector<std::size_t> bulk_element_ids =
            findElementsInMesh(bulk_mesh, element_node_ids_bulk);

        if (bulk_element_ids.empty())
        {
            ERR("No element could be found for the subdomain element {:d}. "
                "Corresponding bulk mesh node ids are:",
                e->getID());
            for (auto const i : element_node_ids_bulk)
            {
                ERR("\t{:d}", i);
            }
            OGS_FATAL(
                "Expect at least one element to be found in the bulk mesh.");
        }

        bulk_element_ids_map[e->getID()] = std::move(bulk_element_ids);
    }

    return bulk_element_ids_map;
}

/// Updates or checks the existing mesh's property with the given values.
void updateOrCheckExistingSubdomainProperty(
    MeshLib::Mesh& mesh, std::string const& property_name,
    std::vector<std::size_t> const& values,
    MeshLib::MeshItemType const mesh_item_type, bool const force_overwrite)
{
    auto& properties = mesh.getProperties();
    if (!properties.existsPropertyVector<std::size_t>(property_name))
    {
        addPropertyToMesh(mesh, property_name, mesh_item_type, 1, values);
        return;
    }

    //
    // Check the existing property against new values.
    //
    auto& original_property =
        *properties.getPropertyVector<std::size_t>(property_name);
    if (std::equal(begin(original_property), end(original_property),
                   begin(values), end(values)))
    {
        INFO(
            "There is already a '{:s}' property present in the subdomain mesh "
            "'{:s}' and it is equal to the newly computed values.",
            property_name, mesh.getName());
        return;
    }

    //
    // Property differs. Notify and update if forced.
    //
    WARN(
        "There is already a '{:s}' property present in the subdomain mesh "
        "'{:s}' and it is not equal to the newly computed values.",
        property_name,
        mesh.getName());

    if (!force_overwrite)
    {
        OGS_FATAL("The force overwrite flag was not specified, exiting.");
    }

    INFO("Overwriting '{:s}' property.", property_name);
    original_property.resize(values.size());
    std::copy(begin(values), end(values), begin(original_property));
}
}  // namespace

namespace MeshGeoToolsLib
{
void identifySubdomainMesh(MeshLib::Mesh& subdomain_mesh,
                           MeshLib::Mesh const& bulk_mesh,
                           MeshNodeSearcher const& mesh_node_searcher,
                           bool const force_overwrite = false)
{
    auto const& bulk_node_ids =
        identifySubdomainMeshNodes(subdomain_mesh, mesh_node_searcher);

    updateOrCheckExistingSubdomainProperty(
        subdomain_mesh, "bulk_node_ids", bulk_node_ids,
        MeshLib::MeshItemType::Node, force_overwrite);

    auto const& bulk_element_ids =
        identifySubdomainMeshElements(subdomain_mesh, bulk_mesh);

    // The bulk_element_ids could be of two types: one element per entry---this
    // is the expected case for the boundary meshes; multiple elements per
    // entry---this happens if the subdomain mesh lies inside the bulk mesh and
    // has lower dimension.
    // First find out the type, then add/check the CellData or FieldData.
    if (all_of(begin(bulk_element_ids), end(bulk_element_ids),
               [](std::vector<std::size_t> const& v) { return v.size() == 1; }))
    {
        // All vectors are of size 1, so the data can be flattened and
        // stored in CellData or compared to existing CellData.
        std::vector<std::size_t> unique_bulk_element_ids;
        unique_bulk_element_ids.reserve(bulk_element_ids.size());
        transform(begin(bulk_element_ids), end(bulk_element_ids),
                  back_inserter(unique_bulk_element_ids),
                  [](std::vector<std::size_t> const& v) { return v[0]; });

        updateOrCheckExistingSubdomainProperty(
            subdomain_mesh, "bulk_element_ids", unique_bulk_element_ids,
            MeshLib::MeshItemType::Cell, force_overwrite);
    }
    else
    {
        // Some of the boundary elements are connected to multiple bulk
        // elements; Store the array in FieldData with additional CellData array
        // for the number of elements, which also provides the offsets.
        std::vector<std::size_t> flat_bulk_element_ids;
        flat_bulk_element_ids.reserve(2 * bulk_element_ids.size());  // Guess.
        std::vector<std::size_t> number_of_bulk_element_ids;
        number_of_bulk_element_ids.reserve(bulk_element_ids.size());

        for (std::vector<std::size_t> const& v : bulk_element_ids)
        {
            number_of_bulk_element_ids.push_back(v.size());
            flat_bulk_element_ids.insert(end(flat_bulk_element_ids), begin(v),
                                         end(v));
        }

        updateOrCheckExistingSubdomainProperty(
            subdomain_mesh, "number_bulk_elements", number_of_bulk_element_ids,
            MeshLib::MeshItemType::Cell, force_overwrite);
        updateOrCheckExistingSubdomainProperty(
            subdomain_mesh, "bulk_element_ids", flat_bulk_element_ids,
            MeshLib::MeshItemType::IntegrationPoint, force_overwrite);
    }
}
}  // namespace MeshGeoToolsLib
