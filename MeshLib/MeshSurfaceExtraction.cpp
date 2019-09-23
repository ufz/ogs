/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshSurfaceExtraction class.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshSurfaceExtraction.h"

#include <memory>

#include <boost/math/constants/constants.hpp>

#include <logog/include/logog.hpp>

#include "MeshLib/Elements/Line.h"
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

    std::vector<MeshLib::Node*> sfc_nodes;
    std::vector<std::size_t> node_id_map(subsfc_mesh.getNumberOfNodes());
    get2DSurfaceNodes(sfc_nodes, subsfc_mesh.getNumberOfNodes(), sfc_elements,
                      node_id_map);

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

    // transmit the original node ids of the subsurface mesh as a property
    if (!subsfc_node_id_prop_name.empty())
    {
        MeshLib::addPropertyToMesh(*result, subsfc_node_id_prop_name,
                                   MeshLib::MeshItemType::Node, 1, id_map);
    }

    // transmit the original subsurface element ids as a property
    if (!subsfc_element_id_prop_name.empty())
    {
        MeshLib::addPropertyToMesh(*result, subsfc_element_id_prop_name,
                                   MeshLib::MeshItemType::Cell, 1,
                                   element_ids_map);
    }

    // transmit the face id of the original subsurface element as a property
    if (!face_id_prop_name.empty())
    {
        MeshLib::addPropertyToMesh(*result, face_id_prop_name,
                                   MeshLib::MeshItemType::Cell, 1,
                                   face_ids_map);
    }

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

void MeshSurfaceExtraction::get2DSurfaceNodes(
    std::vector<MeshLib::Node*>& sfc_nodes, std::size_t n_all_nodes,
    const std::vector<MeshLib::Element*>& sfc_elements,
    std::vector<std::size_t>& node_id_map)
{
    std::vector<const MeshLib::Node*> tmp_nodes(n_all_nodes, nullptr);
    for (auto const* elem : sfc_elements)
    {
        for (unsigned j = 0; j < elem->getNumberOfBaseNodes(); ++j)
        {
            const MeshLib::Node* node(elem->getNode(j));
            tmp_nodes[node->getID()] = node;
        }
    }
    const std::size_t nNodes(tmp_nodes.size());
    for (unsigned i = 0; i < nNodes; ++i)
    {
        if (tmp_nodes[i])
        {
            node_id_map[i] = sfc_nodes.size();
            sfc_nodes.push_back(new MeshLib::Node(*tmp_nodes[i]));
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

    std::vector<MeshLib::Node*> sfc_nodes;
    std::vector<std::size_t> node_id_map(mesh.getNumberOfNodes());
    get2DSurfaceNodes(sfc_nodes, mesh.getNumberOfNodes(), sfc_elements,
                      node_id_map);

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

}  // end namespace MeshLib
