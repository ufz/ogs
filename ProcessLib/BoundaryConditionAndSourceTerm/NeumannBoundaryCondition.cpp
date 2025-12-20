// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NeumannBoundaryCondition.h"

#include "ParameterLib/Utils.h"

namespace ProcessLib
{
NeumannBoundaryConditionConfig parseNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "Neumann");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Neumann__parameter}
    auto const name = config.getConfigParameter<std::string>("parameter");
    DBUG("{}: parameter {:s}", __FUNCTION__, name);

    auto const area_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Neumann__area_parameter}
        config.getConfigParameterOptional<std::string>("area_parameter");
    if (area_parameter_name.has_value())
    {
        DBUG("{}: area parameter name '{:s}'", __FUNCTION__,
             area_parameter_name.value());
    }

    return {name, area_parameter_name};
}

std::unique_ptr<NeumannBoundaryCondition> createNeumannBoundaryCondition(
    NeumannBoundaryConditionConfig const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing Neumann BC.");

    auto const& param = ParameterLib::findParameter<double>(
        config.parameter_name, parameters, 1, &bc_mesh);

    // In case of partitioned mesh the boundary could be empty, i.e. there
    // is no boundary condition.
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

    ParameterLib::Parameter<double> const* integral_measure(nullptr);
    if (global_dim - bc_mesh.getDimension() != 1)
    {
        if (!config.area_parameter_name)
        {
            OGS_FATAL("{}: tag <area_parameter> required in Neumann BC.",
                      __FUNCTION__);
        }
        integral_measure = &ParameterLib::findParameter<double>(
            config.area_parameter_name.value(), parameters, 1, &bc_mesh);
    }

    if (bc_mesh.getDimension() >= global_dim)
    {
        OGS_FATAL(
            "The dimension ({:d}) of the given boundary mesh '{:s}' is not "
            "lower than the bulk dimension ({:d}).",
            bc_mesh.getDimension(), bc_mesh.getName(), global_dim);
    }

    return std::make_unique<NeumannBoundaryCondition>(
        integration_order, shapefunction_order, dof_table, variable_id,
        component_id, global_dim, bc_mesh,
        NeumannBoundaryConditionData{param, integral_measure});
}

}  // namespace ProcessLib
