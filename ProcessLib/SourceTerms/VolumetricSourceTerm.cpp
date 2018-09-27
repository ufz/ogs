/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
    NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
    unsigned const integration_order, unsigned const shapefunction_order,
    int const variable_id, int const component_id,
    Parameter<double> const& volumetric_source_term)
    : SourceTerm(source_term_dof_table),
      _volumetric_source_term(volumetric_source_term)
{
    // check basic data consistency
    if (variable_id >=
        static_cast<int>(source_term_dof_table.getNumberOfVariables()))
    {
        OGS_FATAL(
            "Variable id too high. Actual value: %d, maximum value: %d.",
            variable_id,
            source_term_dof_table.getNumberOfVariables());
    }
    if (component_id >=
        source_term_dof_table.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Component id too high. Actual value: %d, maximum value: %d.",
            component_id,
            source_term_dof_table.getNumberOfVariableComponents(variable_id));
    }

    ProcessLib::createLocalAssemblers<VolumetricSourceTermLocalAssembler>(
        source_term_mesh.getDimension(), source_term_mesh.getElements(),
        source_term_dof_table, shapefunction_order, _local_assemblers,
        source_term_mesh.isAxiallySymmetric(), integration_order,
        _volumetric_source_term);
}

void VolumetricSourceTerm::integrate(const double t,
                                     GlobalVector& b) const
{
    DBUG("Assemble VolumetricSourceTerm.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberOnDereferenced(
        &VolumetricSourceTermLocalAssemblerInterface::integrate,
        _local_assemblers, _source_term_dof_table, t, b);
}

}   // namespace ProcessLib
