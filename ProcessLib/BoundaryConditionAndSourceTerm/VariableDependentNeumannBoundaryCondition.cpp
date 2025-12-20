// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "VariableDependentNeumannBoundaryCondition.h"

#include "ParameterLib/Utils.h"

namespace ProcessLib
{
VariableDependentNeumannConfig parseVariableDependentNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "VariableDependentNeumann");

    auto const constant_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__VariableDependentNeumann__constant_name}
        config.getConfigParameter<std::string>("constant_name");

    auto const coefficient_current_variable_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__VariableDependentNeumann__coefficient_current_variable_name}
        config.getConfigParameter<std::string>(
            "coefficient_current_variable_name");

    auto const coefficient_other_variable_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__VariableDependentNeumann__coefficient_other_variable_name}
        config.getConfigParameter<std::string>(
            "coefficient_other_variable_name");

    auto const coefficient_mixed_variables_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__VariableDependentNeumann__coefficient_mixed_variables_name}
        config.getConfigParameter<std::string>(
            "coefficient_mixed_variables_name");

    return {constant_name, coefficient_current_variable_name,
            coefficient_other_variable_name, coefficient_mixed_variables_name};
}

std::unique_ptr<VariableDependentNeumannBoundaryCondition>
createVariableDependentNeumannBoundaryCondition(
    VariableDependentNeumannConfig const& coefficients,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing VariableDependentNeumann BC.");
    if (dof_table.getNumberOfVariables() != 2)
    {
        OGS_FATAL(
            "VariableDependentNeumann BC only implemented for 2 variable "
            "processes.");
    }
    assert(variable_id == 0 || variable_id == 1);

    if (bc_mesh.getDimension() + 1 != global_dim)
    {
        OGS_FATAL(
            "The dimension ({:d}) of the given boundary mesh '{:s}' is not by "
            "one lower than the bulk dimension ({:d}).",
            bc_mesh.getDimension(), bc_mesh.getName(), global_dim);
    }

    auto const& constant = ParameterLib::findParameter<double>(
        coefficients.constant_name, parameters, 1, &bc_mesh);

    auto const& coefficient_current_variable =
        ParameterLib::findParameter<double>(coefficients.current_variable_name,
                                            parameters, 1, &bc_mesh);

    auto const& coefficient_other_variable =
        ParameterLib::findParameter<double>(coefficients.other_variable_name,
                                            parameters, 1, &bc_mesh);

    auto const& coefficient_mixed_variables =
        ParameterLib::findParameter<double>(coefficients.mixed_variables_name,
                                            parameters, 1, &bc_mesh);

    std::vector<MeshLib::Node*> const& bc_nodes = bc_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(bc_mesh, bc_nodes);
    auto dof_table_boundary_other_variable =
        dof_table.deriveBoundaryConstrainedMap(
            (variable_id + 1) % 2, {component_id}, std::move(bc_mesh_subset));

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

    return std::make_unique<VariableDependentNeumannBoundaryCondition>(
        integration_order, shapefunction_order, dof_table, variable_id,
        component_id, global_dim, bc_mesh,
        VariableDependentNeumannBoundaryConditionData{
            constant, coefficient_current_variable, coefficient_other_variable,
            coefficient_mixed_variables,
            std::move(dof_table_boundary_other_variable)});
}

}  // namespace ProcessLib
