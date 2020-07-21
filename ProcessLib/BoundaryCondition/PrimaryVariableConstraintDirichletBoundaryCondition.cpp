/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PrimaryVariableConstraintDirichletBoundaryCondition.h"

#include <algorithm>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
PrimaryVariableConstraintDirichletBoundaryCondition::
    PrimaryVariableConstraintDirichletBoundaryCondition(
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        ParameterLib::Parameter<double> const& threshold_parameter,
        bool const less)
    : _parameter(parameter),
      _bc_mesh(bc_mesh),
      _variable_id(variable_id),
      _component_id(component_id),
      _threshold_parameter(threshold_parameter),
      _less(less)
{
    checkParametersOfDirichletBoundaryCondition(_bc_mesh, dof_table_bulk,
                                                _variable_id, _component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = bc_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset)));
}

void PrimaryVariableConstraintDirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    ParameterLib::SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    auto const& nodes_in_bc_mesh = _bc_mesh.getNodes();
    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + nodes_in_bc_mesh.size());
    bc_values.values.reserve(bc_values.values.size() + nodes_in_bc_mesh.size());
    for (auto const* const node : nodes_in_bc_mesh)
    {
        auto const id = node->getID();
        // TODO: that might be slow, but only done once
        auto const global_index = _dof_table_boundary->getGlobalIndex(
            {_bc_mesh.getID(), MeshLib::MeshItemType::Node, id}, _variable_id,
            _component_id);
        if (global_index == NumLib::MeshComponentMap::nop)
        {
            continue;
        }
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply the
        // Dirichlet BC, the negative index is not accepted like other matrix or
        // vector PETSc routines. Therefore, the following if-condition is
        // applied.
        if (global_index >= 0)
        {
            // fetch the value of the primary variable
            auto const local_x = x.get(std::vector{global_index});
            pos.setNodeID(id);
            pos.setCoordinates(*node);
            if (_less && local_x[0] < _threshold_parameter(t, pos).front())
            {
                bc_values.ids.emplace_back(global_index);
                bc_values.values.emplace_back(_parameter(t, pos).front());
            }
            else if (!_less &&
                     local_x[0] > _threshold_parameter(t, pos).front())
            {
                bc_values.ids.emplace_back(global_index);
                bc_values.values.emplace_back(_parameter(t, pos).front());
            }
        }
    }
}

std::unique_ptr<PrimaryVariableConstraintDirichletBoundaryCondition>
createPrimaryVariableConstraintDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing PrimaryVariableConstraintDirichletBoundaryCondition from "
        "config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "PrimaryVariableConstraintDirichlet");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__PrimaryVariableConstraintDirichletBoundaryCondition__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter {:s}", param_name);

    auto& parameter = ParameterLib::findParameter<double>(
        param_name, parameters, 1, &bc_mesh);

    auto const threshold_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__PrimaryVariableConstraintDirichletBoundaryCondition__threshold_parameter}
        config.getConfigParameter<std::string>("threshold_parameter");
    DBUG("Using parameter {:s} as threshold_parameter",
         threshold_parameter_name);

    auto& threshold_parameter = ParameterLib::findParameter<double>(
        threshold_parameter_name, parameters, 1, &bc_mesh);

    auto const comparison_operator_string =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__PrimaryVariableConstraintDirichletBoundaryCondition__comparison_operator}
        config.getConfigParameter<std::string>("comparison_operator");
    if (comparison_operator_string != "greater" &&
        comparison_operator_string != "less")
    {
        OGS_FATAL(
            "The comparison operator is '{:s}', but has to be either "
            "'greater' or 'less'.",
            comparison_operator_string);
    }
    bool const less = comparison_operator_string == "less";

// In case of partitioned mesh the boundary could be empty, i.e. there is no
// boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createBoundaryCondition().
    if (bc_mesh.getDimension() == 0 && bc_mesh.getNumberOfNodes() == 0 &&
        bc_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<
        PrimaryVariableConstraintDirichletBoundaryCondition>(
        parameter, bc_mesh, dof_table_bulk, variable_id, component_id,
        threshold_parameter, less);
}

}  // namespace ProcessLib
