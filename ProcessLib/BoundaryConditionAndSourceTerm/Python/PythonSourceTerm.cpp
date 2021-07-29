/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonSourceTerm.h"

#include <pybind11/pybind11.h>

#include <iostream>

#include "FlushStdoutGuard.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Utils/CreateLocalAssemblers.h"
#include "PythonSourceTermLocalAssembler.h"

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{
PythonSourceTerm::PythonSourceTerm(
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    PythonSourceTermData&& source_term_data, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    bool const flush_stdout)
    : SourceTerm(std::move(source_term_dof_table)),
      _source_term_data(std::move(source_term_data)),
      _flush_stdout(flush_stdout)
{
    BoundaryConditionAndSourceTerm::createLocalAssemblers<
        PythonSourceTermLocalAssembler>(
        global_dim, _source_term_data.source_term_mesh.getElements(),
        *_source_term_dof_table, shapefunction_order, _local_assemblers,
        _source_term_data.source_term_mesh.isAxiallySymmetric(),
        integration_order, _source_term_data);
}

void PythonSourceTerm::integrate(const double t, const GlobalVector& x,
                                 GlobalVector& b, GlobalMatrix* Jac) const
{
    FlushStdoutGuard guard(_flush_stdout);

    GlobalExecutor::executeMemberOnDereferenced(
        &PythonSourceTermLocalAssemblerInterface::assemble, _local_assemblers,
        *_source_term_dof_table, t, x, b, Jac);
}

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
