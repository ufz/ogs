/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonuniformVariableDependantNeumannBoundaryCondition.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace
{
MeshLib::PropertyVector<double> const* getAndCheckPropertyVector(
    std::string const& name, MeshLib::Mesh const& mesh)
{
    auto const* const property_vector =
        mesh.getProperties().getPropertyVector<double>(name);
    if (property_vector == nullptr)
    {
        OGS_FATAL("A property with name `%s' does not exist in `%s'.",
                  name.c_str(), mesh.getName().c_str());
    }
    if (property_vector->getMeshItemType() != MeshLib::MeshItemType::Node)
    {
        OGS_FATAL(
            "Only nodal fields are supported for non-uniform fields. Field "
            "`%s' is not nodal.",
            name.c_str());
    }
    if (property_vector->getNumberOfComponents() != 1)
    {
        OGS_FATAL("`%s' is not a one-component field.", name.c_str());
    }
    return property_vector;
};
}  // namespace

namespace ProcessLib
{
std::unique_ptr<NonuniformVariableDependantNeumannBoundaryCondition>
createNonuniformVariableDependantNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, MeshLib::Mesh const& bulk_mesh)
{
    DBUG("Constructing NonuniformVariableDependantNeumann BC from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NonuniformVariableDependantNeumann");
    if (dof_table.getNumberOfVariables() != 2)
    {
        OGS_FATAL(
            "NonuniformVariableDependantNeumann BC only implemented for 2 "
            "variable processes.");
    }
    assert(variable_id == 0 || variable_id == 1);

    // TODO finally use ProcessLib::Parameter here
    auto const constant_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependantNeumann__constant_name}
        config.getConfigParameter<std::string>("constant_name");
    //    auto const* const constant = xxx(constant_name);
    auto const* const constant =
        getAndCheckPropertyVector(constant_name, boundary_mesh);

    auto const prefac_current_variable_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependantNeumann__prefac_current_variable_name}
        config.getConfigParameter<std::string>("prefac_current_variable_name");
    auto const* const prefac_current_variable =
        getAndCheckPropertyVector(prefac_current_variable_name, boundary_mesh);

    auto const prefac_other_variable_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependantNeumann__prefac_other_variable_name}
        config.getConfigParameter<std::string>("prefac_other_variable_name");
    auto const* const prefac_other_variable =
        getAndCheckPropertyVector(prefac_other_variable_name, boundary_mesh);

    auto const prefac_mixed_variables_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformVariableDependantNeumann__prefac_mixed_variables_name}
        config.getConfigParameter<std::string>("prefac_mixed_variables_name");
    auto const* const prefac_mixed_variables =
        getAndCheckPropertyVector(prefac_mixed_variables_name, boundary_mesh);

    std::string const mapping_to_bulk_nodes_property = "bulk_node_ids";
    auto const* const mapping_to_bulk_nodes =
        boundary_mesh.getProperties().getPropertyVector<std::size_t>(
            mapping_to_bulk_nodes_property);

    if (!(mapping_to_bulk_nodes && mapping_to_bulk_nodes->getMeshItemType() ==
                                       MeshLib::MeshItemType::Node) &&
        mapping_to_bulk_nodes->getNumberOfComponents() == 1)
    {
        OGS_FATAL("Field `%s' is not set up properly.",
                  mapping_to_bulk_nodes_property.c_str());
    }

    std::vector<MeshLib::Node*> const& bc_nodes = boundary_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(boundary_mesh, bc_nodes);
    auto const* const dof_table_boundary_other_variable =
        dof_table.deriveBoundaryConstrainedMap(
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
        NonuniformVariableDependantNeumannBoundaryCondition>(
        integration_order, shapefunction_order, dof_table, variable_id,
        component_id, bulk_mesh.getDimension(), boundary_mesh,
        NonuniformVariableDependantNeumannBoundaryConditionData{
            *constant, *prefac_current_variable, *prefac_other_variable,
            *prefac_mixed_variables, *dof_table_boundary_other_variable});
}

}  // namespace ProcessLib
