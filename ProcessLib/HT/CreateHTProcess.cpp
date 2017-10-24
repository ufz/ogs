/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHTProcess.h"

#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"
#include "MaterialLib/PorousMedium/CreatePorousMediaProperties.h"

#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "HTProcess.h"
#include "HTProcessData.h"

namespace ProcessLib
{
namespace HT
{
std::unique_ptr<Process> createHTProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HT");

    DBUG("Create HTProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HT__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {
        //! \ogs_file_param_special{prj__processes__process__HT__process_variables__temperature}
        "temperature",
        //! \ogs_file_param_special{prj__processes__process__HT__process_variables__pressure}
        "pressure"});

    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties{
        MaterialLib::PorousMedium::createPorousMediaProperties(
            mesh, config, parameters)};

    //! \ogs_file_param{prj__processes__process__HT__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");
    auto fluid_properties =
        MaterialLib::Fluid::createFluidProperties(fluid_config);

    // Parameter for the density of the fluid.
    auto& fluid_reference_density= findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__fluid_reference_density}
        "fluid_reference_density", parameters, 1);
    DBUG("Use \'%s\' as fluid_reference_density parameter.",
         fluid_reference_density.name.c_str());

    // Parameter for the density of the solid.
    auto& density_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__density_solid}
        "density_solid", parameters, 1);
    DBUG("Use \'%s\' as density_solid parameter.", density_solid.name.c_str());

    // Parameter for the specific heat capacity of the solid.
    auto& specific_heat_capacity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__specific_heat_capacity_solid}
        "specific_heat_capacity_solid", parameters, 1);
    DBUG("Use \'%s\' as specific_heat_capacity_solid parameter.",
         specific_heat_capacity_solid.name.c_str());

    // Parameter for the thermal conductivity of the solid (only one scalar per
    // element, i.e., the isotropic case is handled at the moment)
    auto& thermal_dispersivity_longitudinal = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__thermal_dispersivity_longitudinal}
        "thermal_dispersivity_longitudinal", parameters, 1);
    DBUG("Use \'%s\' as thermal_dispersivity_longitudinal parameter.",
         thermal_dispersivity_longitudinal.name.c_str());

    // Parameter for the thermal conductivity of the solid (only one scalar per
    // element, i.e., the isotropic case is handled at the moment)
    auto& thermal_dispersivity_transversal = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__thermal_dispersivity_transversal}
        "thermal_dispersivity_transversal", parameters, 1);
    DBUG("Use \'%s\' as thermal_dispersivity_transversal parameter.",
         thermal_dispersivity_transversal.name.c_str());

    // Parameter for the thermal conductivity of the solid (only one scalar per
    // element, i.e., the isotropic case is handled at the moment)
    auto& thermal_conductivity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__thermal_conductivity_solid}
        "thermal_conductivity_solid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_solid parameter.",
         thermal_conductivity_solid.name.c_str());

    // Parameter for the thermal conductivity of the fluid.
    auto& thermal_conductivity_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__thermal_conductivity_fluid}
        "thermal_conductivity_fluid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_fluid parameter.",
         thermal_conductivity_fluid.name.c_str());

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__HT__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(b.size() > 0 && b.size() < 4);
    if (b.size() < mesh.getDimension())
        OGS_FATAL(
            "specific body force (gravity vector) has %d components, mesh "
            "dimension is %d",
            b.size(), mesh.getDimension());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    HTProcessData process_data{std::move(porous_media_properties),
                               density_solid,
                               fluid_reference_density,
                               std::move(fluid_properties),
                               thermal_dispersivity_longitudinal,
                               thermal_dispersivity_transversal,
                               specific_heat_capacity_solid,
                               thermal_conductivity_solid,
                               thermal_conductivity_fluid,
                               specific_body_force,
                               has_gravity};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HT_temperature_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::make_unique<HTProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}

}  // namespace HT
}  // namespace ProcessLib
