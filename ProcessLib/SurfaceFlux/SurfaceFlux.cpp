/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SurfaceFlux.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
SurfaceFlux::SurfaceFlux(
    MeshLib::Mesh& boundary_mesh,
    std::size_t bulk_property_number_of_components,
    unsigned const integration_order)
{
    DBUG("Create local balance assemblers.");
    // Populate the vector of local assemblers.
    _local_assemblers.resize(boundary_mesh.getElements().size());

    // needed to create dof table
    auto mesh_subset_all_nodes = std::make_unique<MeshLib::MeshSubset>(
        boundary_mesh, boundary_mesh.getNodes());

    // Collect the mesh subsets in a vector.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    bulk_property_number_of_components,
                    [&]() { return *mesh_subset_all_nodes; });

    // needed for creation of local assemblers
    auto dof_table = std::make_unique<NumLib::LocalToGlobalIndexMap const>(
        std::move(all_mesh_subsets), NumLib::ComponentOrder::BY_LOCATION);

    auto const bulk_element_ids =
        boundary_mesh.getProperties().template getPropertyVector<std::size_t>(
            "bulk_element_ids");
    auto const bulk_face_ids =
        boundary_mesh.getProperties().template getPropertyVector<std::size_t>(
            "bulk_face_ids");

    ProcessLib::createLocalAssemblers<SurfaceFluxLocalAssembler>(
        boundary_mesh.getDimension() + 1,  // or bulk_mesh.getDimension()?
        boundary_mesh.getElements(), *dof_table, 1, _local_assemblers,
        boundary_mesh.isAxiallySymmetric(), integration_order,
        *bulk_element_ids, *bulk_face_ids);
}

void SurfaceFlux::integrate(
    GlobalVector const& x,
    MeshLib::PropertyVector<double>& balance,
    double const t,
    MeshLib::Mesh const& bulk_mesh,
    std::function<Eigen::Vector3d(std::size_t const, MathLib::Point3d const&,
                                  double const, GlobalVector const&)> const&
        getFlux)
{
    DBUG("Integrate SurfaceFlux.");

    GlobalExecutor::executeMemberOnDereferenced(
        &SurfaceFluxLocalAssemblerInterface::integrate,
        _local_assemblers, x, balance, t, bulk_mesh, getFlux);
}

}  // namespace ProcessLib
