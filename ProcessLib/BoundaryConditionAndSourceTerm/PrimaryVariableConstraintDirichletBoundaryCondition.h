/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
/// The PrimaryVariableConstraintDirichletBoundaryCondition is a Dirichlet-type
/// boundary condition.
///
/// The class implements a constraint Dirichlet-type boundary condition. Using
/// a threshold for the primary variable which is given by a parameter the
/// boundary condition can be switched on/off. The value that is set a
/// Dirichlet-type boundary condition is given by another parameter.
class PrimaryVariableConstraintDirichletBoundaryCondition final
    : public BoundaryCondition
{
public:
    PrimaryVariableConstraintDirichletBoundaryCondition(
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        ParameterLib::Parameter<double> const& threshold_parameter,
        bool const less);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    ///< parameter that defines the Dirirchlet-type condition values
    ParameterLib::Parameter<double> const& _parameter;

    MeshLib::Mesh const& _bc_mesh;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;

    /// The threshold parameter used to the switch on/off the Dirichlet-type
    /// boundary condition.
    ParameterLib::Parameter<double> const& _threshold_parameter;

    /// The value less is used for the calculation of the constraint
    /// criterion. If less is set to true (i.e. 'less' is set in the
    /// project file) the criterion 'calculated_value
    /// < _threshold_parameter' is evaluated to switch on/off the boundary
    /// condition.
    /// If less will be set to false in case 'greater' is given in the project
    /// file and the condition 'calculated_value > _threshold_parameter' is
    /// evaluated.
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
