/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VolumetricSourceTerm.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
VolumetricSourceTerm::VolumetricSourceTerm(
    MeshLib::Mesh const& source_term_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    unsigned const integration_order, unsigned const shapefunction_order,
    ParameterLib::Parameter<double> const& volumetric_source_term)
    : SourceTerm(std::move(source_term_dof_table)),
      _volumetric_source_term(volumetric_source_term)
{
    ProcessLib::createLocalAssemblers<VolumetricSourceTermLocalAssembler>(
        source_term_mesh.getDimension(), source_term_mesh.getElements(),
        *_source_term_dof_table, shapefunction_order, _local_assemblers,
        source_term_mesh.isAxiallySymmetric(), integration_order,
        _volumetric_source_term);
}

void VolumetricSourceTerm::integrate(const double t, GlobalVector const& /*x*/,
                                     GlobalVector& b,
                                     GlobalMatrix* /*jac*/) const
{
    DBUG("Assemble VolumetricSourceTerm.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberOnDereferenced(
        &VolumetricSourceTermLocalAssemblerInterface::integrate,
        _local_assemblers, *_source_term_dof_table, t, b);
}

}   // namespace ProcessLib
