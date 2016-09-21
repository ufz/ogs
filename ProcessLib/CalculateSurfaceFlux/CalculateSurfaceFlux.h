/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_CALCULATESURFACEFLUX_H_
#define PROCESS_LIB_CALCULATESURFACEFLUX_H_

#include "ProcessLib/Process.h"
#include "CalculateSurfaceFluxLocalAssembler.h"

namespace ProcessLib
{
class CalculateSurfaceFlux final
{
public:
    /// Constructs an object that is able to compute a balance.
    /// @param mesh This mesh represents the boundary that is integrated over.
    /// @param bulk_property_number_of_components The number of components the
    /// variable has.
    CalculateSurfaceFlux(MeshLib::Mesh& boundary_mesh,
                   std::size_t bulk_property_number_of_components);

    /// Executes for each element of the mesh the local intergration procedure.
    /// @param x The global solution the intergration values will be fetched of.
    /// @param balance The vector the integration results will be stored in.
    /// @param bulk_process Stores the variable that is used for the
    /// integration.
    void integrate(GlobalVector const& x,
                   MeshLib::PropertyVector<double>& balance,
                   Process const& bulk_process);

private:
    // the local assemblers for each element of the mesh
    std::vector<std::unique_ptr<CalculateSurfaceFluxLocalAssemblerInterface>>
        _local_assemblers;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_CALCULATESURFACEFLUX_H_
