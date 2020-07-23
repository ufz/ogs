/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LineSourceTerm.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
LineSourceTerm::LineSourceTerm(
    unsigned const bulk_mesh_dimension, MeshLib::Mesh const& source_term_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    unsigned const integration_order, unsigned const shapefunction_order,
    ParameterLib::Parameter<double> const& source_term_parameter)
    : SourceTerm(std::move(source_term_dof_table)),
      _source_term_parameter(source_term_parameter)
{
    ProcessLib::createLocalAssemblers<LineSourceTermLocalAssembler>(
        bulk_mesh_dimension, source_term_mesh.getElements(),
        *_source_term_dof_table, shapefunction_order, _local_assemblers,
        source_term_mesh.isAxiallySymmetric(), integration_order,
        _source_term_parameter);
}

void LineSourceTerm::integrate(const double t, GlobalVector const& /*x*/,
                               GlobalVector& b, GlobalMatrix* /*jac*/) const
{
    DBUG("Assemble LineSourceTerm.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberOnDereferenced(
        &LineSourceTermLocalAssemblerInterface::integrate, _local_assemblers,
        *_source_term_dof_table, t, b);
}

}  // namespace ProcessLib
