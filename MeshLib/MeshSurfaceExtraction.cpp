/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshSurfaceExtraction class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshSurfaceExtraction.h"

#include <boost/math/constants/constants.hpp>
#include <logog/include/logog.hpp>
#include <memory>

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

    INFO("Total surface Area: %f", total_area);

    return node_area_vec;
}

MeshLib::Mesh* MeshSurfaceExtraction::getMeshSurface(
    const MeshLib::Mesh& subsfc_mesh, const MathLib::Vector3& dir, double angle,
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

    auto [sfc_nodes, node_id_map] =
        createNodesFromElements(sfc_elements, subsfc_mesh.getNumberOfNodes());

    // create new elements vector with newly created nodes (and delete
    // temp-elements)
    std::vector<MeshLib::Element*> new_elements;
    try
    {
        new_elements =
            createSfcElementVector(sfc_elements, sfc_nodes, node_id_map);
    }
    catch (std::runtime_error const& err)
    {
        ERR("MeshSurfaceExtraction; could not create new surface "
            "elements:\n%s.",
            err.what());
        std::for_each(sfc_elements.begin(), sfc_elements.end(),
                      [](MeshLib::Element* e) { delete e; });
        std::for_each(sfc_nodes.begin(), sfc_nodes.end(),
                      [](MeshLib::Node* n) { delete n; });
        return nullptr;
    }
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

MeshLib::Mesh* MeshSurfaceExtraction::getMeshBoundary(const MeshLib::Mesh& mesh)
{
    if (mesh.getDimension() == 1)
    {
        return nullptr;
    }

    // For 3D meshes return the 2D surface
    if (mesh.getDimension() == 3)
    {
        MathLib::Vector3 dir(0, 0, 0);
        return getMeshSurface(mesh, dir, 90);
    }

    // For 2D meshes return the boundary lines
    std::vector<MeshLib::Node*> nodes =
        MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> boundary_elements;

    std::vector<MeshLib::Element*> const& org_elems(mesh.getElements());
    for (auto elem : org_elems)
    {
        std::size_t const n_edges(elem->getNumberOfEdges());
        for (std::size_t i = 0; i < n_edges; ++i)
        {
            if (elem->getNeighbor(i) == nullptr)
            {
                MeshLib::Element const* const edge(elem->getEdge(i));
                boundary_elements.push_back(MeshLib::copyElement(edge, nodes));
                delete edge;
            }
        }
    }
    MeshLib::Mesh* result =
        new MeshLib::Mesh("Boundary Mesh", nodes, boundary_elements);
    MeshLib::NodeSearch ns(*result);
    if (ns.searchUnused() == 0)
    {
        return result;
    }
    auto removed = MeshLib::removeNodes(*result, ns.getSearchedNodeIDs(),
                                        result->getName());
    delete result;
    return removed;
}

void MeshSurfaceExtraction::get2DSurfaceElements(
    const std::vector<MeshLib::Element*>& all_elements,
    std::vector<MeshLib::Element*>& sfc_elements,
    std::vector<std::size_t>& element_to_bulk_element_id_map,
    std::vector<std::size_t>& element_to_bulk_face_id_map,
    const MathLib::Vector3& dir, double angle, unsigned mesh_dimension)
{
    if (mesh_dimension < 2 || mesh_dimension > 3)
        ERR("Cannot handle meshes of dimension %i", mesh_dimension);

    bool const complete_surface = (MathLib::scalarProduct(dir, dir) == 0);

    double const pi(boost::math::constants::pi<double>());
    double const cos_theta(std::cos(angle * pi / 180.0));
    MathLib::Vector3 const norm_dir(dir.getNormalizedVector());

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
                if (MathLib::scalarProduct(
                        FaceRule::getSurfaceNormal(face).getNormalizedVector(),
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
                    if (MathLib::scalarProduct(
                            FaceRule::getSurfaceNormal(face.get())
                                .getNormalizedVector(),
                            norm_dir) < cos_theta)
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
    const MeshLib::Mesh& mesh, const MathLib::Vector3& dir, double angle)
{
    INFO("Extracting surface nodes...");
    std::vector<MeshLib::Element*> sfc_elements;
    std::vector<std::size_t> element_to_bulk_element_id_map;
    std::vector<std::size_t> element_to_bulk_face_id_map;
    get2DSurfaceElements(
        mesh.getElements(), sfc_elements, element_to_bulk_element_id_map,
        element_to_bulk_face_id_map, dir, angle, mesh.getDimension());

    auto [sfc_nodes, node_id_map] =
        createNodesFromElements(sfc_elements, mesh.getNumberOfNodes());

    for (auto e : sfc_elements)
    {
        delete e;
    }

    return sfc_nodes;
}

std::vector<MeshLib::Element*> MeshSurfaceExtraction::createSfcElementVector(
    std::vector<MeshLib::Element*> const& sfc_elements,
    std::vector<MeshLib::Node*> const& sfc_nodes,
    std::vector<std::size_t> const& node_id_map)
{
    std::vector<MeshLib::Element*> new_elements;
    new_elements.reserve(sfc_elements.size());
    for (auto sfc_element : sfc_elements)
    {
        unsigned const n_elem_nodes(sfc_element->getNumberOfBaseNodes());
        auto** new_nodes = new MeshLib::Node*[n_elem_nodes];
        for (unsigned k(0); k < n_elem_nodes; k++)
        {
            new_nodes[k] =
                sfc_nodes[node_id_map[sfc_element->getNode(k)->getID()]];
        }
        switch (sfc_element->getGeomType())
        {
            case MeshElemType::TRIANGLE:
                new_elements.push_back(new MeshLib::Tri(new_nodes));
                break;
            case MeshElemType::QUAD:
                new_elements.push_back(new MeshLib::Quad(new_nodes));
                break;
            case MeshElemType::LINE:
                new_elements.push_back(new MeshLib::Line(new_nodes));
                break;
            default:
                OGS_FATAL(
                    "MeshSurfaceExtraction::createSfcElementVector Unknown "
                    "element type '%s'.",
                    MeshElemType2String(sfc_element->getGeomType()).c_str());
        }
    }
    return new_elements;
}

bool MeshSurfaceExtraction::createSfcMeshProperties(
    MeshLib::Mesh& sfc_mesh,
    MeshLib::Properties const& properties,
    std::vector<std::size_t> const& node_ids_map,
    std::vector<std::size_t> const& element_ids_map)
{
    std::size_t const n_elems(sfc_mesh.getNumberOfElements());
    std::size_t const n_nodes(sfc_mesh.getNumberOfNodes());
    if (node_ids_map.size() != n_nodes)
    {
        ERR("MeshSurfaceExtraction::createSfcMeshProperties() - Incorrect "
            "number of node IDs (%d) compared to actual number of surface "
            "nodes "
            "(%d).",
            node_ids_map.size(), n_nodes);
        return false;
    }

    if (element_ids_map.size() != n_elems)
    {
        ERR("MeshSurfaceExtraction::createSfcMeshProperties() - Incorrect "
            "number of element IDs (%d) compared to actual number of surface "
            "elements (%d).",
            element_ids_map.size(), n_elems);
        return false;
    }

    std::size_t vectors_copied(0);
    std::size_t vectors_skipped(0);
    std::vector<std::string> const& array_names =
        properties.getPropertyVectorNames();
    for (std::string const& name : array_names)
    {
        if (processPropertyVector<double>(name, MeshLib::MeshItemType::Cell,
                                          properties, n_elems, element_ids_map,
                                          sfc_mesh) ||
            processPropertyVector<int>(name, MeshLib::MeshItemType::Cell,
                                       properties, n_elems, element_ids_map,
                                       sfc_mesh) ||
            processPropertyVector<double>(name, MeshLib::MeshItemType::Node,
                                          properties, n_nodes, node_ids_map,
                                          sfc_mesh) ||
            processPropertyVector<int>(name, MeshLib::MeshItemType::Node,
                                       properties, n_nodes, node_ids_map,
                                       sfc_mesh))
        {
            vectors_copied++;
        }
        else
        {
            WARN("Skipping property vector '%s' - no matching data type found.",
                 name.c_str());
            vectors_skipped++;
        }
    }
    INFO("%d property vectors copied, %d vectors skipped.", vectors_copied,
         vectors_skipped);
    return true;
}

std::tuple<std::vector<MeshLib::Node*>, std::vector<std::size_t>>
MeshSurfaceExtraction::createNodesFromElements(
    std::vector<MeshLib::Element*> const& elements,
    std::size_t const n_all_nodes)
{
    std::vector<const MeshLib::Node*> tmp_nodes(n_all_nodes, nullptr);
    for (auto const* elem : elements)
    {
        for (unsigned j = 0; j < elem->getNumberOfBaseNodes(); ++j)
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

void MeshSurfaceExtraction::createSurfaceElementsFromElement(
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

        auto const face = std::unique_ptr<MeshLib::Element const>{
            surface_element.getBoundary(j)};
        switch (face->getGeomType())
        {
            case MeshElemType::TRIANGLE:
            {
                surface_elements.push_back(new MeshLib::Tri(
                    *static_cast<const MeshLib::Tri*>(face.get())));
                break;
            }
            case MeshElemType::QUAD:
            {
                surface_elements.push_back(new MeshLib::Quad(
                    *static_cast<const MeshLib::Quad*>(face.get())));
                break;
            }
            case MeshElemType::LINE:
            {
                surface_elements.push_back(new MeshLib::Line(
                    *static_cast<const MeshLib::Line*>(face.get())));
                break;
            }
            case MeshElemType::POINT:
            {
                surface_elements.push_back(new MeshLib::Point(
                    *static_cast<const MeshLib::Point*>(face.get())));
                break;
            }
            default:
                ERR("Unknown face element extracted.");
        }
        element_to_bulk_face_id_map.push_back(j);
        element_to_bulk_element_id_map.push_back(surface_element.getID());
    }
}

std::tuple<std::vector<MeshLib::Element*>, std::vector<std::size_t>,
           std::vector<std::size_t>>
MeshSurfaceExtraction::createSurfaceElements(MeshLib::Mesh const& bulk_mesh)
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

std::unique_ptr<MeshLib::Mesh> MeshSurfaceExtraction::getBoundaryElementsAsMesh(
    MeshLib::Mesh const& bulk_mesh,
    std::string const& subsfc_node_id_prop_name,
    std::string const& subsfc_element_id_prop_name,
    std::string const& face_id_prop_name)
{
    auto const mesh_dimension = bulk_mesh.getDimension();
    if (mesh_dimension < 2 || mesh_dimension > 3)
    {
        ERR("Cannot handle meshes of dimension %i", mesh_dimension);
    }

    // create surface elements based on the subsurface nodes
    std::vector<std::size_t> element_to_bulk_element_id_map;
    std::vector<std::size_t> element_to_bulk_face_id_map;
    std::vector<MeshLib::Element*> surface_elements;
    std::tie(surface_elements, element_to_bulk_element_id_map,
             element_to_bulk_face_id_map) = createSurfaceElements(bulk_mesh);

    // create new nodes needed for the new surface elements
    std::vector<MeshLib::Node*> surface_nodes;
    std::vector<std::size_t> node_id_map;
    std::tie(surface_nodes, node_id_map) =
        MeshSurfaceExtraction::createNodesFromElements(
            surface_elements, bulk_mesh.getNumberOfNodes());

    // create new elements using newly created nodes and delete temp-elements
    std::vector<MeshLib::Element*> new_elements;
    try
    {
        new_elements = MeshSurfaceExtraction::createSfcElementVector(
            surface_elements, surface_nodes, node_id_map);
    }
    catch (std::runtime_error const& err)
    {
        ERR("BoundaryExtraction; could not create new surface elements:\n%s.",
            err.what());
        std::for_each(surface_elements.begin(), surface_elements.end(),
                      [](MeshLib::Element* e) { delete e; });
        std::for_each(surface_nodes.begin(), surface_nodes.end(),
                      [](MeshLib::Node* n) { delete n; });
        return nullptr;
    }
    std::for_each(surface_elements.begin(), surface_elements.end(),
                  [](MeshLib::Element* e) { delete e; });

    std::vector<std::size_t> nodes_to_bulk_nodes_id_map;
    nodes_to_bulk_nodes_id_map.reserve(surface_nodes.size());
    std::transform(begin(surface_nodes), end(surface_nodes),
                   std::back_inserter(nodes_to_bulk_nodes_id_map),
                   [](MeshLib::Node* const n) { return n->getID(); });

    std::unique_ptr<MeshLib::Mesh> surface_mesh(new Mesh(
        bulk_mesh.getName() + "-Surface", surface_nodes, new_elements));

    addBulkIDPropertiesToMesh(
        *surface_mesh, subsfc_node_id_prop_name, nodes_to_bulk_nodes_id_map,
        subsfc_element_id_prop_name, element_to_bulk_element_id_map,
        face_id_prop_name, element_to_bulk_face_id_map);

    /// Copies relevant parts of scalar arrays to the surface mesh
    if (!MeshSurfaceExtraction::createSfcMeshProperties(
            *surface_mesh, bulk_mesh.getProperties(),
            nodes_to_bulk_nodes_id_map, element_to_bulk_element_id_map))
    {
        ERR("Transferring subsurface properties failed.");
    }
    return surface_mesh;
}

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
