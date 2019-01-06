/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonuniformVariableDependentNeumannBoundaryCondition.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
std::unique_ptr<NonuniformVariableDependentNeumannBoundaryCondition>
createNonuniformVariableDependentNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, MeshLib::Mesh const& bulk_mesh)
{
    DBUG("Constructing NonuniformVariableDependentNeumann BC from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NonuniformVariableDependentNeumann");
    if (dof_table.getNumberOfVariables() != 2)
    {
        OGS_FATAL(
            "NonuniformVariableDependentNeumann BC only implemented for 2 "
            "variable processes.");
    }
    assert(variable_id == 0 || variable_id == 1);

    // TODO finally use ProcessLib::Parameter here
    auto const constant_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependentNeumann__constant_name}
        config.getConfigParameter<std::string>("constant_name");
    auto const& constant =
        *boundary_mesh.getProperties().getPropertyVector<double>(
            constant_name, MeshLib::MeshItemType::Node, 1);

    auto const coefficient_current_variable_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependentNeumann__coefficient_current_variable_name}
        config.getConfigParameter<std::string>(
            "coefficient_current_variable_name");
    auto const& coefficient_current_variable =
        *boundary_mesh.getProperties().getPropertyVector<double>(
            coefficient_current_variable_name, MeshLib::MeshItemType::Node, 1);

    auto const coefficient_other_variable_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependentNeumann__coefficient_other_variable_name}
        config.getConfigParameter<std::string>(
            "coefficient_other_variable_name");
    auto const& coefficient_other_variable =
        *boundary_mesh.getProperties().getPropertyVector<double>(
            coefficient_other_variable_name, MeshLib::MeshItemType::Node, 1);

    auto const coefficient_mixed_variables_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependentNeumann__coefficient_mixed_variables_name}
        config.getConfigParameter<std::string>(
            "coefficient_mixed_variables_name");
    auto const& coefficient_mixed_variables =
        *boundary_mesh.getProperties().getPropertyVector<double>(
            coefficient_mixed_variables_name, MeshLib::MeshItemType::Node, 1);

    std::string const mapping_to_bulk_nodes_property = "bulk_node_ids";
    boundary_mesh.getProperties().getPropertyVector<std::size_t>(
        mapping_to_bulk_nodes_property, MeshLib::MeshItemType::Node, 1);

    std::vector<MeshLib::Node*> const& bc_nodes = boundary_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(boundary_mesh, bc_nodes);
    auto const& dof_table_boundary_other_variable =
        *dof_table.deriveBoundaryConstrainedMap(
            (variable_id + 1) % 2, {component_id}, std::move(bc_mesh_subset));

    // In case of partitioned mesh the boundary could be empty, i.e. there is no
    // boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createBoundaryCondition().
    if (boundary_mesh.getDimension() == 0 &&
        boundary_mesh.getNumberOfNodes() == 0 &&
        boundary_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<
        NonuniformVariableDependentNeumannBoundaryCondition>(
        integration_order, shapefunction_order, dof_table, variable_id,
        component_id, bulk_mesh.getDimension(), boundary_mesh,
        NonuniformVariableDependentNeumannBoundaryConditionData{
            constant, coefficient_current_variable, coefficient_other_variable,
            coefficient_mixed_variables, dof_table_boundary_other_variable});
}

}  // namespace ProcessLib
