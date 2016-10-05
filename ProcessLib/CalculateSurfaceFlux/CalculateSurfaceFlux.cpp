/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CalculateSurfaceFlux.h"

#include <cassert>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
CalculateSurfaceFlux::CalculateSurfaceFlux(
    MeshLib::Mesh& boundary_mesh,
    std::size_t bulk_property_number_of_components,
    unsigned const integration_order)
{
    DBUG("Create local balance assemblers.");
    // Populate the vector of local assemblers.
    _local_assemblers.resize(boundary_mesh.getElements().size());

    // needed to create dof table
    std::unique_ptr<MeshLib::MeshSubset const> mesh_subset_all_nodes(
        new MeshLib::MeshSubset(boundary_mesh, &boundary_mesh.getNodes()));

    // Collect the mesh subsets in a vector.
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
    std::generate_n(
        std::back_inserter(all_mesh_subsets),
        bulk_property_number_of_components,
        [&]() {
            return std::unique_ptr<MeshLib::MeshSubsets>{
                new MeshLib::MeshSubsets{mesh_subset_all_nodes.get()}};
        });

    // needed for creation of local assemblers
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> dof_table(
        new NumLib::LocalToGlobalIndexMap(std::move(all_mesh_subsets),
                                          NumLib::ComponentOrder::BY_LOCATION));

    auto const bulk_element_ids =
        boundary_mesh.getProperties().template getPropertyVector<std::size_t>(
            "OriginalSubsurfaceElementIDs");
    if (!bulk_element_ids)
        OGS_FATAL(
            "OriginalSubsurfaceElementIDs boundary mesh property not found.");
    auto const bulk_face_ids =
        boundary_mesh.getProperties().template getPropertyVector<std::size_t>(
            "OriginalFaceIDs");
    if (!bulk_face_ids)
        OGS_FATAL("OriginalFaceIDs boundary mesh property not found.");

    ProcessLib::createLocalAssemblers<CalculateSurfaceFluxLocalAssembler>(
        boundary_mesh.getDimension() + 1,  // or bulk_mesh.getDimension()?
        boundary_mesh.getElements(), *dof_table, _local_assemblers,
        boundary_mesh.isAxiallySymmetric(), integration_order,
        *bulk_element_ids, *bulk_face_ids);
}

void CalculateSurfaceFlux::integrate(GlobalVector const& x,
                               MeshLib::PropertyVector<double>& balance,
                               Process const& bulk_process)
{
    DBUG("Integrate CalculateSurfaceFlux.");

    GlobalExecutor::executeMemberOnDereferenced(
        &CalculateSurfaceFluxLocalAssemblerInterface::integrate, _local_assemblers, x,
        balance, bulk_process);
}

}  // namespace ProcessLib
