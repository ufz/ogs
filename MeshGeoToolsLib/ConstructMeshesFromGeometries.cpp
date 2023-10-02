/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ConstructMeshesFromGeometries.h"

#ifdef USE_PETSC
#include <mpi.h>

#include "MeshLib/NodePartitionedMesh.h"
#include "MeshLib/Utils/transformMeshToNodePartitionedMesh.h"
#endif

#include "BaseLib/Logging.h"
#include "BoundaryElementsSearcher.h"
#include "GeoLib/GEOObjects.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"
#include "MeshLib/Utils/createMeshFromElementSelection.h"
#include "MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{
std::string meshNameFromGeometry(std::string const& geometrical_set_name,
                                 std::string const& geometry_name)
{
    return geometrical_set_name + "_" + geometry_name;
}

template <typename GeometryVec>
std::vector<std::unique_ptr<MeshLib::Mesh>>
constructAdditionalMeshesFromGeometries(
    std::vector<GeometryVec*> const& geometries,
    MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher,
    bool const multiple_nodes_allowed)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> additional_meshes;

    for (GeometryVec* const geometry_vec : geometries)
    {
        // Each geometry_vec has a name, this is the first part of the full
        // name.
        auto const& vec_name = geometry_vec->getName();

        auto const& vec_data = geometry_vec->getVector();

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

            DBUG("Creating mesh from geometry {:s} {:s}.", vec_name,
                 geometry_name);

            auto subdomain_mesh = createMeshFromElementSelection(
                meshNameFromGeometry(vec_name, geometry_name),
                MeshLib::cloneElements(
                    boundary_element_searcher.getBoundaryElements(
                        geometry, multiple_nodes_allowed)));

#ifdef USE_PETSC
            // The subdomain_mesh is not yet a NodePartitionedMesh.
            // The bulk_mesh, which is a NodePartitionedMesh, is needed to
            // construct the subdomain NodePartitionedMesh
            auto const* bulk_mesh =
                dynamic_cast<MeshLib::NodePartitionedMesh const*>(
                    &boundary_element_searcher.mesh);

            additional_meshes.push_back(
                MeshLib::transformMeshToNodePartitionedMesh(
                    bulk_mesh, subdomain_mesh.get()));
#else
            // Nothing special to be done in the serial case.
            additional_meshes.emplace_back(std::move(subdomain_mesh));
#endif
        }
    }
    return additional_meshes;
}

std::vector<std::unique_ptr<MeshLib::Mesh>>
constructAdditionalMeshesFromGeoObjects(GeoLib::GEOObjects const& geo_objects,
                                        MeshLib::Mesh const& mesh,
                                        std::unique_ptr<SearchLength>
                                            search_length_algorithm,
                                        bool const multiple_nodes_allowed)
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
            geo_objects.getPoints(), boundary_element_searcher,
            multiple_nodes_allowed);
        std::move(begin(point_meshes), end(point_meshes),
                  std::back_inserter(additional_meshes));
    }

    //
    // Polylines
    //
    {
        auto polyline_meshes = constructAdditionalMeshesFromGeometries(
            geo_objects.getPolylines(), boundary_element_searcher,
            multiple_nodes_allowed);
        std::move(begin(polyline_meshes), end(polyline_meshes),
                  std::back_inserter(additional_meshes));
    }

    // Surfaces
    {
        auto surface_meshes = constructAdditionalMeshesFromGeometries(
            geo_objects.getSurfaces(), boundary_element_searcher,
            multiple_nodes_allowed);
        std::move(begin(surface_meshes), end(surface_meshes),
                  std::back_inserter(additional_meshes));
    }

    // Set axial symmetry for the additional meshes to the same value as the
    // "bulk" mesh.
    std::for_each(begin(additional_meshes), end(additional_meshes),
                  [axial_symmetry = mesh.isAxiallySymmetric()](auto& m)
                  { m->setAxiallySymmetric(axial_symmetry); });
    return additional_meshes;
}
}  // namespace MeshGeoToolsLib
