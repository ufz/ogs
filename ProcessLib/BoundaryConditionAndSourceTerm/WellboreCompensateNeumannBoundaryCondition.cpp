/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "WellboreCompensateNeumannBoundaryCondition.h"

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
WellboreCompensateCoefficients parseWellboreCompensateNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config)
{
    DBUG("Parsing WellboreCompensateNeumann BC.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "WellboreCompensateNeumann");

    auto const pressure_coefficient =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__WellboreCompensateNeumann__coefficient_pressure}
        config.getConfigParameter<double>("coefficient_pressure");

    auto const velocity_coefficient =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__WellboreCompensateNeumann__coefficient_velocity}
        config.getConfigParameter<double>("coefficient_velocity");

    auto const enthalpy_coefficient =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__WellboreCompensateNeumann__coefficient_enthalpy}
        config.getConfigParameter<double>("coefficient_enthalpy");

    return {pressure_coefficient, velocity_coefficient, enthalpy_coefficient};
}

std::unique_ptr<WellboreCompensateNeumannBoundaryCondition>
createWellboreCompensateNeumannBoundaryCondition(
    WellboreCompensateCoefficients const& coefficients,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    DBUG("Constructing WellboreCompensateNeumann BC.");

    if (dof_table.getNumberOfVariables() != 3)
    {
        OGS_FATAL(
            "WellboreCompensateNeumann BC only implemented for the "
            "WellboreSimulator processes.");
    }
    assert(variable_id == 0 || variable_id == 1 || variable_id == 2);

    if (bc_mesh.getDimension() + 1 != global_dim)
    {
        OGS_FATAL(
            "The dimension ({:d}) of the given boundary mesh '{:s}' is not by "
            "one "
            "lower than the bulk dimension ({:d}).",
            bc_mesh.getDimension(), bc_mesh.getName(), global_dim);
    }

    auto const pressure_id = 0;
    auto const velocity_id = 1;
    auto const enthalpy_id = 2;

    std::vector<MeshLib::Node*> const& bc_nodes = bc_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(bc_mesh, bc_nodes);
    auto dof_table_boundary_pressure = dof_table.deriveBoundaryConstrainedMap(
        pressure_id, {component_id}, std::move(bc_mesh_subset));
    auto dof_table_boundary_velocity = dof_table.deriveBoundaryConstrainedMap(
        velocity_id, {component_id}, std::move(bc_mesh_subset));
    auto dof_table_boundary_enthalpy = dof_table.deriveBoundaryConstrainedMap(
        enthalpy_id, {component_id}, std::move(bc_mesh_subset));

    //
    // maybe the boundary mesh needs material ids
    auto media_map = MaterialPropertyLib::createMaterialSpatialDistributionMap(
        media, bc_mesh);

    return std::make_unique<WellboreCompensateNeumannBoundaryCondition>(
        integration_order, shapefunction_order, dof_table, variable_id,
        component_id, global_dim, bc_mesh,
        WellboreCompensateNeumannBoundaryConditionData{
            coefficients, std::move(dof_table_boundary_pressure),
            std::move(dof_table_boundary_velocity),
            std::move(dof_table_boundary_enthalpy), std::move(media_map)});
}

}  // namespace ProcessLib
