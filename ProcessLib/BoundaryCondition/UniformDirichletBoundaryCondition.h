/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_BOUNDARY_CONDITION_H_
#define PROCESS_LIB_BOUNDARY_CONDITION_H_

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "BoundaryCondition.h"

namespace ProcessLib
{
/// The UniformDirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class UniformDirichletBoundaryCondition : public BoundaryCondition
{
public:
    UniformDirichletBoundaryCondition(
        NumLib::IndexValueVector<GlobalIndexType>&& bc)
        : _bc(std::move(bc))
    {
    }

    void apply(const double, GlobalVector const&, GlobalMatrix&,
               GlobalVector&) override
    {
    }

    void getDirichletBCValues(
        const double /*t*/,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
    {
        if (!_bc.ids.empty())
            bc_values = std::move(_bc);
    }

private:
    mutable NumLib::IndexValueVector<GlobalIndexType> _bc;
};

std::unique_ptr<UniformDirichletBoundaryCondition>
createUniformDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, std::vector<std::size_t>&& mesh_node_ids,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t const mesh_id,
    int const variable_id, int const component_id);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
