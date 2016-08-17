/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_DIRICHLETBOUNDARYCONDITION_H
#define PROCESSLIB_DIRICHLETBOUNDARYCONDITION_H

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "BoundaryCondition.h"

namespace ProcessLib
{
/// The DirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class DirichletBoundaryCondition : public BoundaryCondition
{
public:
    DirichletBoundaryCondition(
        NumLib::IndexValueVector<GlobalIndexType>&& bc)
        : _bc(std::move(bc))
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

std::unique_ptr<DirichletBoundaryCondition>
createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, std::vector<std::size_t>&& mesh_node_ids,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t const mesh_id,
    int const variable_id, int const component_id);

}  // namespace ProcessLib

#endif  // PROCESSLIB_DIRICHLETBOUNDARYCONDITION_H
