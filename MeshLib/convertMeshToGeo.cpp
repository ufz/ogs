/**
 * \file
 * \author Karsten Rink
 * \date   2013-07-05
 * \brief  Implementation of  of mesh to geometry conversion.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "convertMeshToGeo.h"

#include <logog/include/logog.hpp>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Triangle.h"
#include "GeoLib/Surface.h"

#include "Mesh.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "MeshInformation.h"
#include "MeshEditing/MeshRevision.h"

namespace MeshLib {

bool convertMeshToGeo(const MeshLib::Mesh &mesh, GeoLib::GEOObjects &geo_objects, double eps)
{
    if (mesh.getDimension() != 2)
    {
        ERR ("Mesh to geometry conversion is only working for 2D meshes.");
        return false;
    }

    // nodes to points conversion
    std::string mesh_name(mesh.getName());
    {
        auto points = std::make_unique<std::vector<GeoLib::Point*>>();
        points->reserve(mesh.getNumberOfNodes());

        for (auto node_ptr : mesh.getNodes())
            points->push_back(new GeoLib::Point(*node_ptr, node_ptr->getID()));

        geo_objects.addPointVec(std::move(points), mesh_name, nullptr, eps);
    }
    const std::vector<std::size_t> id_map (geo_objects.getPointVecObj(mesh_name)->getIDMap());

    // elements to surface triangles conversion
    std::string const mat_name ("MaterialIDs");
    auto bounds (MeshInformation::getValueBounds<int>(mesh, mat_name));
    const unsigned nMatGroups(bounds.second-bounds.first+1);
    auto sfcs = std::make_unique<std::vector<GeoLib::Surface*>>();
    sfcs->reserve(nMatGroups);
    auto const& points = *geo_objects.getPointVec(mesh_name);
    for (unsigned i=0; i<nMatGroups; ++i)
        sfcs->push_back(new GeoLib::Surface(points));

    const std::vector<MeshLib::Element*> &elements = mesh.getElements();
    const std::size_t nElems (mesh.getNumberOfElements());

    MeshLib::PropertyVector<int> const*const materialIds =
        mesh.getProperties().existsPropertyVector<int>("MaterialIDs")
            ? mesh.getProperties().getPropertyVector<int>("MaterialIDs")
            : nullptr;

    for (unsigned i=0; i<nElems; ++i)
    {
        auto surfaceId = !materialIds ? 0 : ((*materialIds)[i] - bounds.first);
        MeshLib::Element* e (elements[i]);
        if (e->getGeomType() == MeshElemType::TRIANGLE)
            (*sfcs)[surfaceId]->addTriangle(id_map[e->getNodeIndex(0)], id_map[e->getNodeIndex(1)], id_map[e->getNodeIndex(2)]);
        if (e->getGeomType() == MeshElemType::QUAD)
        {
            (*sfcs)[surfaceId]->addTriangle(id_map[e->getNodeIndex(0)], id_map[e->getNodeIndex(1)], id_map[e->getNodeIndex(2)]);
            (*sfcs)[surfaceId]->addTriangle(id_map[e->getNodeIndex(0)], id_map[e->getNodeIndex(2)], id_map[e->getNodeIndex(3)]);
        }
        // all other element types are ignored (i.e. lines)
    }

    std::for_each(sfcs->begin(), sfcs->end(), [](GeoLib::Surface* sfc){ if (sfc->getNumberOfTriangles()==0) delete sfc; sfc = nullptr;});
    auto sfcs_end = std::remove(sfcs->begin(), sfcs->end(), nullptr);
    sfcs->erase(sfcs_end, sfcs->end());

    geo_objects.addSurfaceVec(std::move(sfcs), mesh_name);
    return true;
}

MeshLib::Mesh* convertSurfaceToMesh(const GeoLib::Surface &sfc, const std::string &mesh_name, double eps)
{
    // convert to a mesh including duplicated nodes
    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;
    std::size_t nodeId = 0;
    for (std::size_t i=0; i<sfc.getNumberOfTriangles(); i++)
    {
        auto* tri = sfc[i];
        auto** tri_nodes = new MeshLib::Node*[3];
        for (unsigned j=0; j<3; j++)
            tri_nodes[j] = new MeshLib::Node(tri->getPoint(j)->getCoords(), nodeId++);
        elements.push_back(new MeshLib::Tri(tri_nodes, i));
        for (unsigned j=0; j<3; j++)
            nodes.push_back(tri_nodes[j]);
    }
    MeshLib::Mesh mesh_with_duplicated_nodes(mesh_name, nodes, elements);

    // remove duplicated nodes
    MeshLib::MeshRevision rev(mesh_with_duplicated_nodes);
    return rev.simplifyMesh(mesh_with_duplicated_nodes.getName(), eps);
}

}

