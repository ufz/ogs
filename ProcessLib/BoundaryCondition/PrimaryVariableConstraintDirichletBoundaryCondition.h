/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BoundaryCondition.h"

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib
{
// TODO docu
/// The PrimaryVariableConstraintDirichletBoundaryCondition class describes a
/// constant in space and time PrimaryVariableConstraintDirichlet boundary
/// condition. The expected parameter in the passed configuration is "value"
/// which, when not present defaults to zero.
class PrimaryVariableConstraintDirichletBoundaryCondition final
    : public BoundaryCondition
{
public:
    PrimaryVariableConstraintDirichletBoundaryCondition(
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        ParameterLib::Parameter<double> const& constraint_threshold_parameter,
        bool const less);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    ParameterLib::Parameter<double> const& _parameter;

    MeshLib::Mesh const& _bc_mesh;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;

    /// The threshold value used to the switch off/on the Dirichlet-type
    /// boundary condition.
    ParameterLib::Parameter<double> const& _constraint_threshold_parameter;

    /// The boolean value lower is used for the calculation of the constraint
    /// criterion, i.e., if lower is set to true the criterion 'calculated_value
    /// < constraint_threshold' is evaluated to switch on/off the boundary
    /// condition, else 'calculated_value > constraint_threshold' is evaluated.
    bool const _less;
};

std::unique_ptr<PrimaryVariableConstraintDirichletBoundaryCondition>
createPrimaryVariableConstraintDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>&
        parameters);

}  // namespace ProcessLib
