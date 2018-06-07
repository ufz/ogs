/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Process.h"
#include "CalculateSurfaceFluxLocalAssembler.h"

namespace ProcessLib
{
class CalculateSurfaceFlux final
{
public:
    /// Constructs an object that is able to compute a balance.
    /// @param boundary_mesh The integration boundary mesh.
    /// @param bulk_property_number_of_components The number of components the
    /// variable has.
    /// @param integration_order Integration order used in local assembly.
    CalculateSurfaceFlux(MeshLib::Mesh& boundary_mesh,
                         std::size_t bulk_property_number_of_components,
                         unsigned const integration_order);

    /// Executes for each element of the mesh the local intergration procedure.
    /// @param x The global solution the intergration values will be fetched of.
    /// @param balance The vector the integration results will be stored in.
    /// @param t The balance will be computed at the time t.
    /// @param bulk_mesh Stores a reference to the bulk mesh that is needed to
    /// fetch the information for the integration over the surface mesh.
    /// @param getFlux function that calculates the flux in the integration
    /// points of the face elements of the bulk element that belongs to the
    /// surface.
    void integrate(GlobalVector const& x,
                   MeshLib::PropertyVector<double>& balance,
                   double const t,
                   MeshLib::Mesh const& bulk_mesh,
                   std::function<Eigen::Vector3d(
                       std::size_t const, MathLib::Point3d const&, double const,
                       GlobalVector const&)> const& getFlux);

private:
    // the local assemblers for each element of the mesh
    std::vector<std::unique_ptr<CalculateSurfaceFluxLocalAssemblerInterface>>
        _local_assemblers;
};

}   // namespace ProcessLib
