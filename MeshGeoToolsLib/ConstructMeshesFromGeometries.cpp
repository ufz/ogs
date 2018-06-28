/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ConstructMeshesFromGeometries.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Node.h"

#include "BoundaryElementsSearcher.h"
#include "MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{
template <typename GeometryVec>
std::vector<std::unique_ptr<MeshLib::Mesh>>
constructAdditionalMeshesFromGeometries(
    std::vector<GeometryVec*> const& geometries,
    MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> additional_meshes;

    for (GeometryVec* const geometry_vec : geometries)
    {
        // Each geometry_vec has a name, this is the first part of the full
        // name.
        auto const& vec_name = geometry_vec->getName();

        auto const& vec_data = *geometry_vec->getVector();

        auto const vec_size = geometry_vec->size();
        for (std::size_t i = 0; i < vec_size; ++i)
        {
            // Each geometry has a name, this is the second part of the full
            // name.
            std::string geometry_name;
            bool const is_geometry_named =
                geometry_vec->getNameOfElementByID(i, geometry_name);
            if (!is_geometry_named)
            {
                continue;
            }

            auto const& geometry = *vec_data[i];

            DBUG("Creating mesh from geometry %s %s.", vec_name.c_str(),
                 geometry_name.c_str());

            additional_meshes.emplace_back(createMeshFromElementSelection(
                vec_name + "_" + geometry_name,
                MeshLib::cloneElements(
                    boundary_element_searcher.getBoundaryElements(geometry))));
        }
    }
    return additional_meshes;
}

std::vector<std::unique_ptr<MeshLib::Mesh>>
constructAdditionalMeshesFromGeoObjects(GeoLib::GEOObjects const& geo_objects,
                                        MeshLib::Mesh const& mesh,
                                        std::unique_ptr<SearchLength>
                                            search_length_algorithm)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> additional_meshes;

    auto const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh, std::move(search_length_algorithm));

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    //
    // Points
    //
    {
        auto point_meshes = constructAdditionalMeshesFromGeometries(
            geo_objects.getPoints(), boundary_element_searcher);
        std::move(begin(point_meshes), end(point_meshes),
                  std::back_inserter(additional_meshes));
    }

    //
    // Polylines
    //
    {
        auto polyline_meshes = constructAdditionalMeshesFromGeometries(
            geo_objects.getPolylines(), boundary_element_searcher);
        std::move(begin(polyline_meshes), end(polyline_meshes),
                  std::back_inserter(additional_meshes));
    }

    // Surfaces
    {
        auto surface_meshes = constructAdditionalMeshesFromGeometries(
            geo_objects.getSurfaces(), boundary_element_searcher);
        std::move(begin(surface_meshes), end(surface_meshes),
                  std::back_inserter(additional_meshes));
    }

    // Set axial symmetry for the additional meshes to the same value as the
    // "bulk" mesh.
    std::for_each(begin(additional_meshes), end(additional_meshes),
                  [axial_symmetry = mesh.isAxiallySymmetric()](auto& m) {
                      m->setAxiallySymmetric(axial_symmetry);
                  });
    return additional_meshes;
}
}  // namespace MeshGeoToolsLib
