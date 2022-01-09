/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoRichardsFlowProcess.h"

#include <cassert>

#include "CreateSimplifiedElasticityModel.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "SimplifiedElasticityModel.h"
#include "ThermoRichardsFlowProcess.h"
#include "ThermoRichardsFlowProcessData.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::permeability, MaterialPropertyLib::porosity,
        MaterialPropertyLib::biot_coefficient,
        MaterialPropertyLib::relative_permeability,
        MaterialPropertyLib::saturation};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::viscosity,
        MaterialPropertyLib::density,
    };
    std::array const required_solid_properties = {MaterialPropertyLib::density};

    // Thermal properties are not checked because they can be phase property or
    // meduim property (will be enabled later).
    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, required_medium_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
        checkRequiredProperties(m.second->phase("Solid"),
                                required_solid_properties);
    }
}

void checkProcessVariableComponents(ProcessVariable const& variable)
{
    if (variable.getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from one: got {:d}.",
            variable.getName(),
            variable.getNumberOfGlobalComponents());
    }
}

std::unique_ptr<Process> createThermoRichardsFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_RICHARDS_FLOW");
    DBUG("Create ThermoRichardsFlowProcess.");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

    /// \section processvariablestrf Process Variables

    //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_T;
    ProcessVariable* variable_p;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        /// Primary process variables as they appear in the global component
        /// vector:
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__THERMO_RICHARDS_FLOW__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__THERMO_RICHARDS_FLOW__process_variables__pressure}
             "pressure"});
        variable_T = &per_process_variables[0].get();
        variable_p = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL(
            "So far, only the monolithic scheme is implemented for "
            "THERMO_RICHARDS_FLOW");
    }

    checkProcessVariableComponents(*variable_T);
    checkProcessVariableComponents(*variable_p);

    /// \section parameterstrf Process Parameters
    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != mesh.getDimension())
        {
            OGS_FATAL(
                "specific body force (gravity vector) has {:d} components, "
                "but mesh dimension is {:d}",
                b.size(), mesh.getDimension());
        }
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    DBUG(
        "Check the media properties of ThermoRichardsFlow process "
        "...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    bool const mass_lumping =
        //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__mass_lumping}
        config.getConfigParameter<bool>("mass_lumping", false);

    std::unique_ptr<SimplifiedElasticityModel> simplified_elasticity =
        createElasticityModel(config);

    ThermoRichardsFlowProcessData process_data{
        std::move(media_map), std::move(specific_body_force), mass_lumping,
        std::move(simplified_elasticity)};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<ThermoRichardsFlowProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}
}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
