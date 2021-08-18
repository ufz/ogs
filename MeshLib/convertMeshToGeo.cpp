/**
 * \file
 * \author Karsten Rink
 * \date   2013-07-05
 * \brief  Implementation of  of mesh to geometry conversion.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "convertMeshToGeo.h"

#include "BaseLib/Logging.h"
#include "Elements/Quad.h"
#include "Elements/Tri.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"
#include "Mesh.h"
#include "MeshEditing/MeshRevision.h"
#include "MeshInformation.h"

namespace
{
/// Convert and add mesh nodes to the geo_objects. A new name of the geo object
/// is returned.
std::string convertMeshNodesToGeoPoints(MeshLib::Mesh const& mesh,
                                        double const eps,
                                        GeoLib::GEOObjects& geo_objects)
{
    auto points = std::make_unique<std::vector<GeoLib::Point*>>();
    points->reserve(mesh.getNumberOfNodes());

    std::transform(begin(mesh.getNodes()), end(mesh.getNodes()),
                   std::back_inserter(*points),
                   [](MeshLib::Node const* node_ptr) {
                       return new GeoLib::Point(*node_ptr, node_ptr->getID());
                   });

    auto geoobject_name = mesh.getName();  // Possibly modified when adding
                                           // points to the geo objects.
    geo_objects.addPointVec(std::move(points), geoobject_name, nullptr, eps);
    return geoobject_name;
}

void addElementToSurface(MeshLib::Element const& e,
                         std::vector<std::size_t> const& id_map,
                         GeoLib::Surface& surface)
{
    if (e.getGeomType() == MeshLib::MeshElemType::TRIANGLE)
    {
        surface.addTriangle(id_map[getNodeIndex(e, 0)],
                            id_map[getNodeIndex(e, 1)],
                            id_map[getNodeIndex(e, 2)]);
        return;
    }
    if (e.getGeomType() == MeshLib::MeshElemType::QUAD)
    {
        surface.addTriangle(id_map[getNodeIndex(e, 0)],
                            id_map[getNodeIndex(e, 1)],
                            id_map[getNodeIndex(e, 2)]);
        surface.addTriangle(id_map[getNodeIndex(e, 0)],
                            id_map[getNodeIndex(e, 2)],
                            id_map[getNodeIndex(e, 3)]);
        return;
    }
    // all other element types are ignored (i.e. lines)
};

}  // namespace

namespace MeshLib
{
bool convertMeshToGeo(const MeshLib::Mesh& mesh,
                      GeoLib::GEOObjects& geo_objects,
                      double const eps)
{
    if (mesh.getDimension() != 2)
    {
        ERR("Mesh to geometry conversion is only working for 2D meshes.");
        return false;
    }

    // Special handling of the bounds in case there are no materialIDs present.
    auto get_material_ids_and_bounds =
        [&]() -> std::tuple<MeshLib::PropertyVector<int> const*,
                            std::pair<int, int>> {
        auto const materialIds = materialIDs(mesh);
        if (!materialIds)
        {
            return std::make_tuple(nullptr, std::make_pair(0, 0));
        }

        auto const bounds = MeshInformation::getValueBounds(*materialIds);
        if (!bounds)
        {
            OGS_FATAL(
                "Could not get minimum/maximum ranges values for the "
                "MaterialIDs property in the mesh '{:s}'.",
                mesh.getName());
        }
        return std::make_tuple(materialIds, *bounds);
    };

    auto const [materialIds, bounds] = get_material_ids_and_bounds();
    // elements to surface triangles conversion
    const unsigned nMatGroups(bounds.second - bounds.first + 1);
    auto sfcs = std::make_unique<std::vector<GeoLib::Surface*>>();
    sfcs->reserve(nMatGroups);
    std::string const geoobject_name =
        convertMeshNodesToGeoPoints(mesh, eps, geo_objects);
    auto const& points = *geo_objects.getPointVec(geoobject_name);
    for (unsigned i = 0; i < nMatGroups; ++i)
    {
        sfcs->push_back(new GeoLib::Surface(points));
    }

    const std::vector<std::size_t>& id_map(
        geo_objects.getPointVecObj(geoobject_name)->getIDMap());
    const std::vector<MeshLib::Element*>& elements = mesh.getElements();
    const std::size_t nElems(mesh.getNumberOfElements());

    for (unsigned i = 0; i < nElems; ++i)
    {
        auto surfaceId = !materialIds ? 0 : ((*materialIds)[i] - bounds.first);
        addElementToSurface(*elements[i], id_map, *(*sfcs)[surfaceId]);
    }

    std::for_each(sfcs->begin(), sfcs->end(), [](GeoLib::Surface*& sfc) {
        if (sfc->getNumberOfTriangles() == 0)
        {
            delete sfc;
            sfc = nullptr;
        }
    });
    auto sfcs_end = std::remove(sfcs->begin(), sfcs->end(), nullptr);
    sfcs->erase(sfcs_end, sfcs->end());

    geo_objects.addSurfaceVec(std::move(sfcs), geoobject_name);
    return true;
}

MeshLib::Mesh* convertSurfaceToMesh(const GeoLib::Surface& sfc,
                                    const std::string& mesh_name,
                                    double eps)
{
    // convert to a mesh including duplicated nodes
    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;
    std::size_t nodeId = 0;
    for (std::size_t i = 0; i < sfc.getNumberOfTriangles(); i++)
    {
        auto* tri = sfc[i];
        auto** tri_nodes = new MeshLib::Node*[3];
        for (unsigned j = 0; j < 3; j++)
        {
            tri_nodes[j] =
                new MeshLib::Node(tri->getPoint(j)->getCoords(), nodeId++);
        }
        elements.push_back(new MeshLib::Tri(tri_nodes, i));
        for (unsigned j = 0; j < 3; j++)
        {
            nodes.push_back(tri_nodes[j]);
        }
    }
    MeshLib::Mesh mesh_with_duplicated_nodes(mesh_name, nodes, elements);

    // remove duplicated nodes
    MeshLib::MeshRevision rev(mesh_with_duplicated_nodes);
    return rev.simplifyMesh(mesh_with_duplicated_nodes.getName(), eps);
}

}  // namespace MeshLib
