/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GlobalVectorFromNamedFunction.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
GlobalVectorFromNamedFunction::GlobalVectorFromNamedFunction(
    NumLib::SpecificFunctionCaller&& function_caller,
    MeshLib::Mesh const& mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_single,
    SecondaryVariableContext& context)
    : _function_caller(std::move(function_caller))
    , _mesh(mesh)
    , _dof_table_single(dof_table_single)
    , _context(context)
{
    assert(dof_table_single.getNumberOfComponents() == 1);
}

GlobalVector const& GlobalVectorFromNamedFunction::call(
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result)
{
    result = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {_dof_table_single.dofSizeWithoutGhosts(),
         _dof_table_single.dofSizeWithoutGhosts(),
         &_dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType nnodes = _mesh.getNumberOfNodes();

    auto const n_args = _function_caller.getNumberOfUnboundArguments();
    assert(dof_table.getNumberOfComponents() == n_args);
    std::vector<double> args(n_args);

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id) {
        // TODO maybe fill args via callback mechanism or remove this class entirely.
        // Caution: The order of args will be the same as the order of the
        // components in the global vector!
        for (std::size_t i = 0; i < n_args; ++i) {
            args[i] = NumLib::getNodalValue(x, _mesh, dof_table, node_id, i);
        }

        _context.index = node_id;
        auto const value = _function_caller.call(args);

        // TODO Problems with PETSc? (local vs. global index)
        result->set(node_id, value);
    }

    return *result;
}

}  // namespace ProcessLib
