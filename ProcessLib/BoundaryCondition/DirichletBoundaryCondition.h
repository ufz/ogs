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
#include "ProcessLib/Parameter/Parameter.h"
#include "BoundaryCondition.h"

namespace ProcessLib
{
// TODO docu
/// The DirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class DirichletBoundaryCondition final : public BoundaryCondition
{
public:
    DirichletBoundaryCondition(Parameter<double> const& parameter,
                               std::vector<std::size_t>&& mesh_node_ids,
                               NumLib::LocalToGlobalIndexMap const& dof_table,
                               std::size_t const mesh_id, int const variable_id,
                               int const component_id)
        : _parameter(parameter),
          _mesh_node_ids(std::move(mesh_node_ids)),
          _dof_table(dof_table),
          _mesh_id(mesh_id),
          _variable_id(variable_id),
          _component_id(component_id)
    {
    }

    void getDirichletBCValues(
        const double t,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    Parameter<double> const& _parameter;

    std::vector<std::size_t> const _mesh_node_ids;
    NumLib::LocalToGlobalIndexMap const& _dof_table;
    std::size_t const _mesh_id;
    int const _variable_id;
    int const _component_id;
};

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, std::vector<std::size_t>&& mesh_node_ids,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t const mesh_id,
    int const variable_id, int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib

#endif  // PROCESSLIB_DIRICHLETBOUNDARYCONDITION_H
