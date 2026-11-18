// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "SurfaceFluxLocalAssembler.h"

namespace ProcessLib
{
class SurfaceFlux final
{
public:
    /// Constructs an object that is able to compute a balance.
    /// @param boundary_mesh The integration boundary mesh.
    /// @param bulk_dof_table The d.o.f. table of the bulk process; the
    /// d.o.f. table for the boundary mesh is derived from it.
    /// @param integration_order Integration order used in local assembly.
    SurfaceFlux(MeshLib::Mesh& boundary_mesh,
                NumLib::LocalToGlobalIndexMap const& bulk_dof_table,
                unsigned const integration_order);

    /// Executes for each element of the mesh the local integration procedure.
    /// @param x The global solution the integration values will be fetched of.
    /// @param balance The vector the integration results will be stored in.
    /// @param t The balance will be computed at the time t.
    /// @param bulk_mesh Stores a reference to the bulk mesh that is needed to
    /// fetch the information for the integration over the surface mesh.
    /// @param active_element_ids The IDs of active elements. It is empty if
    ///                           there is no deactivated subdomain.
    /// @param getFlux function that calculates the flux in the integration
    /// points of the face elements of the bulk element that belongs to the
    /// surface.
    void integrate(std::vector<GlobalVector*> const& x,
                   MeshLib::PropertyVector<double>& balance,
                   double const t,
                   MeshLib::Mesh const& bulk_mesh,
                   std::vector<std::size_t> const& active_element_ids,
                   std::function<Eigen::Vector3d(
                       std::size_t const, MathLib::Point3d const&, double const,
                       std::vector<GlobalVector*> const&)> const& getFlux);

private:
    // the local assemblers for each element of the mesh
    std::vector<std::unique_ptr<SurfaceFluxLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace ProcessLib
