/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshSurfaceExtraction class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshSurfaceExtraction.h"

#include <boost/math/constants/constants.hpp>
#include <memory>

#include "BaseLib/Logging.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Point.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshSearch/NodeSearch.h"

namespace MeshLib
{
template <typename T>
void processPropertyVector(MeshLib::PropertyVector<T> const& property,
                           std::vector<std::size_t> const& id_map,
                           MeshLib::Mesh& sfc_mesh)
{
    auto const number_of_components = property.getNumberOfGlobalComponents();

    auto sfc_prop = getOrCreateMeshProperty<T>(
        sfc_mesh, property.getPropertyName(), property.getMeshItemType(),
        number_of_components);
    sfc_prop->clear();
    sfc_prop->reserve(id_map.size());

    for (auto bulk_id : id_map)
    {
        std::copy_n(&property.getComponent(bulk_id, 0 /*component_id*/),
                    number_of_components, back_inserter(*sfc_prop));
    }
}

bool createSfcMeshProperties(MeshLib::Mesh& sfc_mesh,
                             MeshLib::Properties const& properties,
                             std::vector<std::size_t> const& node_ids_map,
                             std::vector<std::size_t> const& element_ids_map)
{
    std::size_t const n_elems(sfc_mesh.getNumberOfElements());
    std::size_t const n_nodes(sfc_mesh.getNumberOfNodes());
    if (node_ids_map.size() != n_nodes)
    {
        ERR("createSfcMeshProperties() - Incorrect number of node IDs ({:d}) "
            "compared to actual number of surface nodes ({:d}).",
            node_ids_map.size(), n_nodes);
        return false;
    }

    if (element_ids_map.size() != n_elems)
    {
        ERR("createSfcMeshProperties() - Incorrect number of element IDs "
            "({:d}) compared to actual number of surface elements ({:d}).",
            element_ids_map.size(), n_elems);
        return false;
    }
    std::map<MeshLib::MeshItemType, std::vector<std::size_t> const*> const
        id_maps = {{MeshItemType::Cell, &element_ids_map},
                   {MeshItemType::Node, &node_ids_map}};

    std::size_t vectors_copied(0);
    std::size_t vectors_skipped(0);
    for (auto [name, property] : properties)
    {
        if (property->getMeshItemType() != MeshItemType::Cell &&
            property->getMeshItemType() != MeshItemType::Node)
        {
            WARN(
                "Skipping property vector '{:s}' - not defined on cells or "
                "nodes.",
                name);
            vectors_skipped++;
            continue;
        }

        auto const& id_map = *id_maps.at(property->getMeshItemType());
        if (auto p = dynamic_cast<PropertyVector<double>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<float>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<int>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<unsigned>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<long>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<long long>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p =
                     dynamic_cast<PropertyVector<unsigned long>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<unsigned long long>*>(
                     property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<std::size_t>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else if (auto p = dynamic_cast<PropertyVector<char>*>(property))
        {
            processPropertyVector(*p, id_map, sfc_mesh);
            vectors_copied++;
        }
        else
        {
            WARN(
                "Skipping property vector '{:s}' - no matching data type "
                "'{:s}' found.",
                name, typeid(*property).name());
            vectors_skipped++;
        }
    }
    INFO("{:d} property vectors copied, {:d} vectors skipped.", vectors_copied,
         vectors_skipped);
    return true;
}

std::tuple<std::vector<MeshLib::Node*>, std::vector<std::size_t>>
createNodesAndIDMapFromElements(std::vector<MeshLib::Element*> const& elements,
                                std::size_t const n_all_nodes)
{
    std::vector<MeshLib::Node const*> tmp_nodes(n_all_nodes, nullptr);
    for (auto const* elem : elements)
    {
        auto const n_nodes = elem->getNumberOfNodes();
        for (unsigned j = 0; j < n_nodes; ++j)
        {
            const MeshLib::Node* node(elem->getNode(j));
            tmp_nodes[node->getID()] = node;
        }
    }

    std::vector<MeshLib::Node*> nodes;
    std::vector<std::size_t> node_id_map(n_all_nodes);
    for (unsigned i = 0; i < n_all_nodes; ++i)
    {
        if (tmp_nodes[i])
        {
            node_id_map[i] = nodes.size();
            nodes.push_back(new MeshLib::Node(*tmp_nodes[i]));
        }
    }
    return {nodes, node_id_map};
}

std::vector<double> MeshSurfaceExtraction::getSurfaceAreaForNodes(
    const MeshLib::Mesh& mesh)
{
    std::vector<double> node_area_vec;
    if (mesh.getDimension() != 2)
    {
        ERR("Error in MeshSurfaceExtraction::getSurfaceAreaForNodes() - Given "
            "mesh is no surface mesh (dimension != 2).");
        return node_area_vec;
    }

    double total_area(0);

    // for each node, a vector containing all the element idget every element
    const std::vector<MeshLib::Node*>& nodes = mesh.getNodes();
    const std::size_t nNodes(mesh.getNumberOfNodes());
    for (std::size_t n = 0; n < nNodes; ++n)
    {
        double node_area(0);

        std::vector<MeshLib::Element*> conn_elems = nodes[n]->getElements();
        const std::size_t nConnElems(conn_elems.size());

        for (std::size_t i = 0; i < nConnElems; ++i)
        {
            const MeshLib::Element* elem(conn_elems[i]);
            const unsigned nElemParts =
                (elem->getGeomType() == MeshElemType::TRIANGLE) ? 3 : 4;
            const double area = conn_elems[i]->getContent() / nElemParts;
            node_area += area;
            total_area += area;
        }

        node_area_vec.push_back(node_area);
    }

    INFO("Total surface Area: {:f}", total_area);

    return node_area_vec;
}

MeshLib::Mesh* MeshSurfaceExtraction::getMeshSurface(
    const MeshLib::Mesh& subsfc_mesh, Eigen::Vector3d const& dir, double angle,
    std::string const& subsfc_node_id_prop_name,
    std::string const& subsfc_element_id_prop_name,
    std::string const& face_id_prop_name)
{
    // allow slightly greater angles than 90 degrees for numerical errors
    if (angle < 0 || angle > 91)
    {
        ERR("Supported angle between 0 and 90 degrees only.");
        return nullptr;
    }

    INFO("Extracting mesh surface...");
    std::vector<MeshLib::Element*> sfc_elements;
    std::vector<std::size_t> element_ids_map;
    std::vector<std::size_t> face_ids_map;
    get2DSurfaceElements(subsfc_mesh.getElements(), sfc_elements,
                         element_ids_map, face_ids_map, dir, angle,
                         subsfc_mesh.getDimension());

    if (sfc_elements.empty())
    {
        return nullptr;
    }

    auto [sfc_nodes, node_id_map] = createNodesAndIDMapFromElements(
        sfc_elements, subsfc_mesh.getNumberOfNodes());

    // create new elements vector with newly created nodes (and delete
    // temp-elements)
    auto new_elements =
        MeshLib::copyElementVector(sfc_elements, sfc_nodes, &node_id_map);
    std::for_each(sfc_elements.begin(), sfc_elements.end(),
                  [](MeshLib::Element* e) { delete e; });

    std::vector<std::size_t> id_map;
    id_map.reserve(sfc_nodes.size());
    std::transform(begin(sfc_nodes), end(sfc_nodes), std::back_inserter(id_map),
                   [](MeshLib::Node* const n) { return n->getID(); });

    MeshLib::Mesh* result(
        new Mesh(subsfc_mesh.getName() + "-Surface", sfc_nodes, new_elements));

    addBulkIDPropertiesToMesh(*result, subsfc_node_id_prop_name, id_map,
                              subsfc_element_id_prop_name, element_ids_map,
                              face_id_prop_name, face_ids_map);

    if (!createSfcMeshProperties(*result, subsfc_mesh.getProperties(), id_map,
                                 element_ids_map))
    {
        ERR("Transferring subsurface properties failed.");
    }
    return result;
}

void MeshSurfaceExtraction::get2DSurfaceElements(
    const std::vector<MeshLib::Element*>& all_elements,
    std::vector<MeshLib::Element*>& sfc_elements,
    std::vector<std::size_t>& element_to_bulk_element_id_map,
    std::vector<std::size_t>& element_to_bulk_face_id_map,
    Eigen::Vector3d const& dir, double angle, unsigned mesh_dimension)
{
    if (mesh_dimension < 2 || mesh_dimension > 3)
    {
        ERR("Cannot handle meshes of dimension {:i}", mesh_dimension);
    }

    bool const complete_surface = (dir.dot(dir) == 0);

    double const pi(boost::math::constants::pi<double>());
    double const cos_theta(std::cos(angle * pi / 180.0));
    Eigen::Vector3d const norm_dir(dir.normalized());

    for (auto const* elem : all_elements)
    {
        const unsigned element_dimension(elem->getDimension());
        if (element_dimension < mesh_dimension)
        {
            continue;
        }

        if (element_dimension == 2)
        {
            if (!complete_surface)
            {
                auto const* face = elem;
                if (FaceRule::getSurfaceNormal(face).normalized().dot(
                        norm_dir) > cos_theta)
                {
                    continue;
                }
            }
            sfc_elements.push_back(elem->clone());
            element_to_bulk_element_id_map.push_back(elem->getID());
            element_to_bulk_face_id_map.push_back(0);
        }
        else
        {
            if (!elem->isBoundaryElement())
            {
                continue;
            }
            const unsigned nFaces(elem->getNumberOfFaces());
            for (unsigned j = 0; j < nFaces; ++j)
            {
                if (elem->getNeighbor(j) != nullptr)
                {
                    continue;
                }

                auto const face =
                    std::unique_ptr<MeshLib::Element const>{elem->getFace(j)};
                if (!complete_surface)
                {
                    if (FaceRule::getSurfaceNormal(face.get())
                            .normalized()
                            .dot(norm_dir) < cos_theta)
                    {
                        continue;
                    }
                }
                if (face->getGeomType() == MeshElemType::TRIANGLE)
                {
                    sfc_elements.push_back(new MeshLib::Tri(
                        *static_cast<const MeshLib::Tri*>(face.get())));
                }
                else
                {
                    sfc_elements.push_back(new MeshLib::Quad(
                        *static_cast<const MeshLib::Quad*>(face.get())));
                }
                element_to_bulk_element_id_map.push_back(elem->getID());
                element_to_bulk_face_id_map.push_back(j);
            }
        }
    }
}

std::vector<MeshLib::Node*> MeshSurfaceExtraction::getSurfaceNodes(
    const MeshLib::Mesh& mesh, Eigen::Vector3d const& dir, double angle)
{
    INFO("Extracting surface nodes...");
    std::vector<MeshLib::Element*> sfc_elements;
    std::vector<std::size_t> element_to_bulk_element_id_map;
    std::vector<std::size_t> element_to_bulk_face_id_map;

    get2DSurfaceElements(
        mesh.getElements(), sfc_elements, element_to_bulk_element_id_map,
        element_to_bulk_face_id_map, dir, angle, mesh.getDimension());

    std::vector<MeshLib::Node*> surface_nodes;
    std::tie(surface_nodes, std::ignore) =
        createNodesAndIDMapFromElements(sfc_elements, mesh.getNumberOfNodes());

    for (auto e : sfc_elements)
    {
        delete e;
    }

    return surface_nodes;
}

void createSurfaceElementsFromElement(
    MeshLib::Element const& surface_element,
    std::vector<MeshLib::Element*>& surface_elements,
    std::vector<std::size_t>& element_to_bulk_element_id_map,
    std::vector<std::size_t>& element_to_bulk_face_id_map)
{
    const unsigned n_faces(surface_element.getNumberOfBoundaries());
    for (unsigned j = 0; j < n_faces; ++j)
    {
        if (surface_element.getNeighbor(j) != nullptr)
        {
            continue;
        }

        surface_elements.push_back(
            const_cast<MeshLib::Element*>(surface_element.getBoundary(j)));
        element_to_bulk_face_id_map.push_back(j);
        element_to_bulk_element_id_map.push_back(surface_element.getID());
    }
}

std::tuple<std::vector<MeshLib::Element*>, std::vector<std::size_t>,
           std::vector<std::size_t>>
createBoundaryElements(MeshLib::Mesh const& bulk_mesh)
{
    std::vector<std::size_t> element_to_bulk_element_id_map;
    std::vector<std::size_t> element_to_bulk_face_id_map;
    std::vector<MeshLib::Element*> surface_elements;

    auto const& bulk_elements = bulk_mesh.getElements();
    auto const mesh_dimension = bulk_mesh.getDimension();

    for (auto const* elem : bulk_elements)
    {
        const unsigned element_dimension(elem->getDimension());
        if (element_dimension < mesh_dimension)
        {
            continue;
        }

        if (!elem->isBoundaryElement())
        {
            continue;
        }
        createSurfaceElementsFromElement(*elem, surface_elements,
                                         element_to_bulk_element_id_map,
                                         element_to_bulk_face_id_map);
    }
    return {surface_elements, element_to_bulk_element_id_map,
            element_to_bulk_face_id_map};
}

namespace BoundaryExtraction
{
std::unique_ptr<MeshLib::Mesh> getBoundaryElementsAsMesh(
    MeshLib::Mesh const& bulk_mesh,
    std::string const& subsfc_node_id_prop_name,
    std::string const& subsfc_element_id_prop_name,
    std::string const& face_id_prop_name)
{
    auto const mesh_dimension = bulk_mesh.getDimension();
    if (mesh_dimension < 2 || mesh_dimension > 3)
    {
        ERR("Cannot handle meshes of dimension {:i}", mesh_dimension);
    }

    // create boundary elements based on the subsurface nodes
    auto [boundary_elements, element_to_bulk_element_id_map,
          element_to_bulk_face_id_map] = createBoundaryElements(bulk_mesh);

    // create new nodes needed for the new surface elements
    auto [boundary_nodes, node_id_map] = createNodesAndIDMapFromElements(
        boundary_elements, bulk_mesh.getNumberOfNodes());

    // create new elements using newly created nodes and delete temp-elements
    auto new_elements = MeshLib::copyElementVector(
        boundary_elements, boundary_nodes, &node_id_map);
    for (auto* e : boundary_elements)
    {
        delete e;
    }

    std::vector<std::size_t> nodes_to_bulk_nodes_id_map;
    nodes_to_bulk_nodes_id_map.reserve(boundary_nodes.size());
    std::transform(begin(boundary_nodes), end(boundary_nodes),
                   std::back_inserter(nodes_to_bulk_nodes_id_map),
                   [](MeshLib::Node* const n) { return n->getID(); });

    std::unique_ptr<MeshLib::Mesh> boundary_mesh(new Mesh(
        bulk_mesh.getName() + "-boundary", boundary_nodes, new_elements));

    addBulkIDPropertiesToMesh(
        *boundary_mesh, subsfc_node_id_prop_name, nodes_to_bulk_nodes_id_map,
        subsfc_element_id_prop_name, element_to_bulk_element_id_map,
        face_id_prop_name, element_to_bulk_face_id_map);

    /// Copies relevant parts of scalar arrays to the surface mesh
    if (!createSfcMeshProperties(*boundary_mesh, bulk_mesh.getProperties(),
                                 nodes_to_bulk_nodes_id_map,
                                 element_to_bulk_element_id_map))
    {
        ERR("Transferring subsurface properties failed.");
    }
    return boundary_mesh;
}

}  // namespace BoundaryExtraction

void addBulkIDPropertiesToMesh(
    MeshLib::Mesh& surface_mesh,
    std::string const& node_to_bulk_node_id_map_name,
    std::vector<std::size_t> const& node_to_bulk_node_id_map,
    std::string const& element_to_bulk_element_id_map_name,
    std::vector<std::size_t> const& element_to_bulk_element_id_map,
    std::string const& element_to_bulk_face_id_map_name,
    std::vector<std::size_t> const& element_to_bulk_face_id_map)
{
    // transmit the original node ids of the bulk mesh as a property
    if (!node_to_bulk_node_id_map_name.empty())
    {
        MeshLib::addPropertyToMesh(surface_mesh, node_to_bulk_node_id_map_name,
                                   MeshLib::MeshItemType::Node, 1,
                                   node_to_bulk_node_id_map);
    }

    // transmit the original bulk element ids as a property
    if (!element_to_bulk_element_id_map_name.empty())
    {
        MeshLib::addPropertyToMesh(
            surface_mesh, element_to_bulk_element_id_map_name,
            MeshLib::MeshItemType::Cell, 1, element_to_bulk_element_id_map);
    }

    // transmit the face id of the original bulk element as a property
    if (!element_to_bulk_face_id_map_name.empty())
    {
        MeshLib::addPropertyToMesh(
            surface_mesh, element_to_bulk_face_id_map_name,
            MeshLib::MeshItemType::Cell, 1, element_to_bulk_face_id_map);
    }
}

}  // end namespace MeshLib
