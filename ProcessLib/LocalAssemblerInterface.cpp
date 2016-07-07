/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalAssemblerInterface.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
void LocalAssemblerInterface::assemble(
    const std::size_t mesh_item_id,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    assembleConcrete(t, local_x, r_c_indices, M, K, b);
}

void LocalAssemblerInterface::assembleJacobian(
    const std::size_t mesh_item_id,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x, GlobalMatrix& Jac)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    assembleJacobianConcrete(t, local_x, r_c_indices, Jac);
}

void LocalAssemblerInterface::assembleJacobianConcrete(
        double const /*t*/, std::vector<double> const& /*local_x*/,
        NumLib::LocalToGlobalIndexMap::RowColumnIndices const& /*indices*/,
        GlobalMatrix& /*Jac*/)
{
    OGS_FATAL(
        "The assembleJacobian() function is not implemented in the local "
        "assembler.");
}

}  // namespace ProcessLib
