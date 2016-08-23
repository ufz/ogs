/**
 * \file   MeshSurfaceExtraction.cpp
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshSurfaceExtraction class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshSurfaceExtraction.h"

#include <boost/math/constants/constants.hpp>

#include <logog/include/logog.hpp>

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"

namespace MeshLib {

std::vector<double> MeshSurfaceExtraction::getSurfaceAreaForNodes(const MeshLib::Mesh &mesh)
{
    std::vector<double> node_area_vec;
    if (mesh.getDimension() != 2)
    {
        ERR ("Error in MeshSurfaceExtraction::getSurfaceAreaForNodes() - Given mesh is no surface mesh (dimension != 2).");
        return node_area_vec;
    }

    double total_area (0);

    // for each node, a vector containing all the element idget every element
    const std::vector<MeshLib::Node*> &nodes = mesh.getNodes();
    const std::size_t nNodes ( mesh.getNumberOfNodes() );
    for (std::size_t n=0; n<nNodes; ++n)
    {
        double node_area (0);

        std::vector<MeshLib::Element*> conn_elems = nodes[n]->getElements();
        const std::size_t nConnElems (conn_elems.size());

        for (std::size_t i=0; i<nConnElems; ++i)
        {
            const MeshLib::Element* elem (conn_elems[i]);
            const unsigned nElemParts = (elem->getGeomType() == MeshElemType::TRIANGLE) ? 3 : 4;
            const double area = conn_elems[i]->getContent() / nElemParts;
            node_area += area;
            total_area += area;
        }

        node_area_vec.push_back(node_area);
    }

    INFO ("Total surface Area: %f", total_area);

    return node_area_vec;
}

MeshLib::Mesh* MeshSurfaceExtraction::getMeshSurface(
    const MeshLib::Mesh& mesh, const MathLib::Vector3& dir, double angle,
    std::string const& subsfc_node_id_backup_prop_name)
{
    // allow slightly greater angles than 90 degrees for numerical errors
    if (angle < 0 || angle > 91)
    {
        ERR ("Supported angle between 0 and 90 degrees only.");
        return nullptr;
    }

    INFO ("Extracting mesh surface...");
    std::vector<MeshLib::Element*> sfc_elements;
    get2DSurfaceElements(mesh.getElements(), sfc_elements, dir, angle, mesh.getDimension());

    if (sfc_elements.empty())
        return nullptr;

    std::vector<MeshLib::Node*> sfc_nodes;
    std::vector<std::size_t> node_id_map(mesh.getNumberOfNodes());
    get2DSurfaceNodes(sfc_nodes, mesh.getNumberOfNodes(), sfc_elements, node_id_map);

    // create new elements vector with newly created nodes
    std::vector<MeshLib::Element*> new_elements;
    new_elements.reserve(sfc_elements.size());
    for (auto elem = sfc_elements.cbegin(); elem != sfc_elements.cend(); ++elem)
    {
        unsigned const n_elem_nodes ((*elem)->getNumberOfBaseNodes());
        MeshLib::Node** new_nodes = new MeshLib::Node*[n_elem_nodes];
        for (unsigned k(0); k<n_elem_nodes; k++)
            new_nodes[k] = sfc_nodes[node_id_map[(*elem)->getNode(k)->getID()]];
        if ((*elem)->getGeomType() == MeshElemType::TRIANGLE)
            new_elements.push_back(new MeshLib::Tri(new_nodes));
        else {
            assert((*elem)->getGeomType() == MeshElemType::QUAD);
            new_elements.push_back(new MeshLib::Quad(new_nodes));
        }
        delete *elem;
    }

    std::vector<int> id_map;
    if (!subsfc_node_id_backup_prop_name.empty())
    {
        id_map.reserve(sfc_nodes.size());
        for (auto node = sfc_nodes.cbegin(); node != sfc_nodes.cend(); ++node)
            id_map.push_back((*node)->getID());
    }
    MeshLib::Mesh* result (new Mesh(mesh.getName()+"-Surface", sfc_nodes, new_elements));
    // transmit the original node ids of the subsurface mesh as a property
    if (!subsfc_node_id_backup_prop_name.empty()) {
        boost::optional<MeshLib::PropertyVector<int>&> orig_node_ids(
            result->getProperties().createNewPropertyVector<int>(
                subsfc_node_id_backup_prop_name , MeshLib::MeshItemType::Node, 1));
        if (orig_node_ids) {
            orig_node_ids->resize(id_map.size());
            std::copy(id_map.cbegin(), id_map.cend(), orig_node_ids->begin());
        }
    }
    return result;
}

MeshLib::Mesh* MeshSurfaceExtraction::getMeshBoundary(const MeshLib::Mesh &mesh)
{
    if (mesh.getDimension()==1)
        return nullptr;

    // For 3D meshes return the 2D surface
    if (mesh.getDimension()==3)
    {
        MathLib::Vector3 dir(0,0,0);
        return getMeshSurface(mesh, dir, 90);
    }

    // For 2D meshes return the boundary lines
    std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> boundary_elements;

    std::vector<MeshLib::Element*> const& org_elems (mesh.getElements());
    for (auto it=org_elems.begin(); it!=org_elems.end(); ++it)
    {
        MeshLib::Element* elem (*it);
        std::size_t const n_edges (elem->getNumberOfEdges());
        for (std::size_t i=0; i<n_edges; ++i)
            if (elem->getNeighbor(i) == nullptr)
            {
                MeshLib::Element const*const edge (elem->getEdge(i));
                boundary_elements.push_back(MeshLib::copyElement(edge, nodes));
                delete edge;
            }
    }
    MeshLib::Mesh* result = new MeshLib::Mesh("Boundary Mesh", nodes, boundary_elements);
    MeshLib::NodeSearch ns(*result);
    if (ns.searchUnused() == 0) {
        return result;
    } else {
        auto removed = MeshLib::removeNodes(*result, ns.getSearchedNodeIDs(), result->getName());
        delete result;
        return removed;
    }
}

void MeshSurfaceExtraction::get2DSurfaceElements(const std::vector<MeshLib::Element*> &all_elements, std::vector<MeshLib::Element*> &sfc_elements, const MathLib::Vector3 &dir, double angle, unsigned mesh_dimension)
{
    if (mesh_dimension<2 || mesh_dimension>3)
        ERR("Cannot handle meshes of dimension %i", mesh_dimension);

    bool const complete_surface = (MathLib::scalarProduct(dir, dir) == 0);

    double const pi (boost::math::constants::pi<double>());
    double const cos_theta (std::cos(angle * pi / 180.0));
    MathLib::Vector3 const norm_dir (dir.getNormalizedVector());

    for (auto elem = all_elements.cbegin(); elem != all_elements.cend(); ++elem)
    {
        const unsigned element_dimension ((*elem)->getDimension());
        if (element_dimension < mesh_dimension)
            continue;

        if (element_dimension == 2)
        {
            if (!complete_surface)
            {
                MeshLib::Element* face = *elem;
                if (MathLib::scalarProduct(FaceRule::getSurfaceNormal(face).getNormalizedVector(), norm_dir) > cos_theta)
                    continue;
            }
            sfc_elements.push_back(*elem);
        }
        else
        {
            if (!(*elem)->isBoundaryElement())
                continue;
            const unsigned nFaces ((*elem)->getNumberOfFaces());
            for (unsigned j=0; j<nFaces; ++j)
            {
                if ((*elem)->getNeighbor(j) != nullptr)
                    continue;

                auto const face = std::unique_ptr<MeshLib::Element const>{(*elem)->getFace(j)};
                if (!complete_surface)
                {
                    if (MathLib::scalarProduct(FaceRule::getSurfaceNormal(face.get()).getNormalizedVector(), norm_dir) < cos_theta)
                    {
                        continue;
                    }
                }
                if (face->getGeomType() == MeshElemType::TRIANGLE)
                    sfc_elements.push_back(new MeshLib::Tri(*static_cast<const MeshLib::Tri*>(face.get())));
                else
                    sfc_elements.push_back(new MeshLib::Quad(*static_cast<const MeshLib::Quad*>(face.get())));
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
    INFO ("Extracting surface nodes...");
    std::vector<MeshLib::Element*> sfc_elements;
    get2DSurfaceElements(mesh.getElements(), sfc_elements, dir, angle,
                         mesh.getDimension());

    std::vector<MeshLib::Node*> sfc_nodes;
    std::vector<std::size_t> node_id_map(mesh.getNumberOfNodes());
    get2DSurfaceNodes(sfc_nodes, mesh.getNumberOfNodes(), sfc_elements,
                      node_id_map);

    for (auto e : sfc_elements)
        delete e;

    return sfc_nodes;
}

} // end namespace MeshLib
