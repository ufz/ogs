// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateThermoHydroMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelationIce.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "NumLib/NumericalStability/CreateNumericalStabilization.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Common/HydraulicProcess/checkVolumeBalanceEquationSetting.h"
#include "ProcessLib/Common/HydroMechanics/CreateInitialStress.h"
#include "ProcessLib/Common/LiquidThermalExpansion/CheckLiquidThermalExpansivity.h"
#include "ProcessLib/Common/ThermoOsmosis/CheckThermoOsmosisProperties.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "ThermoHydroMechanicsProcess.h"
#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <int DisplacementDim>
std::unique_ptr<Process> createThermoHydroMechanicsProcess(
    std::string const& name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_HYDRO_MECHANICS");
    DBUG("Create ThermoHydroMechanicsProcess.");

    auto coupling_scheme = Common::HydroMechanics::parseCouplingScheme(config);

    auto const* const staggered_scheme =
        std::get_if<Common::HydroMechanics::Staggered>(&coupling_scheme);
    bool const staggered = staggered_scheme != nullptr;
    if (staggered_scheme && staggered_scheme->fixed_stress_over_time_step)
    {
        OGS_FATAL(
            "The coupling scheme option fixed_stress_over_time_step is not "
            "supported by the ThermoHydroMechanics process, which fixes the "
            "volume stress rate over the coupling iteration only.");
    }

    /// \section processvariablesthm Process Variables

    //! \ogs_file_param{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    // Process ids of the individual equations. The monolithic scheme keeps
    // all equations in the single process 0. The staggered scheme solves one
    // equation per process, numbered in the order of the process variable
    // lists assembled below: temperature, pressure, displacement.
    int const thermal_process_id = 0;
    int const hydraulic_process_id = staggered ? 1 : 0;
    int const mechanical_process_id = staggered ? 2 : 0;

    ProcessVariable* variable_T;
    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (!staggered)  // monolithic scheme.
    {
        /// Primary process variables as they appear in the global component
        /// vector:
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_T = &per_process_variables[0].get();
        variable_p = &per_process_variables[1].get();
        variable_u = &per_process_variables[2].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name :
             {"temperature"s, "pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_T = &process_variables[0][0].get();
        variable_p = &process_variables[1][0].get();
        variable_u = &process_variables[2][0].get();
    }

    if (variable_T->getShapeFunctionOrder() != 1)
    {
        OGS_FATAL(
            "The shape function order of temperature must be 1 but its input "
            "value in <process_variable><order> is {:d}. Please correct it.",
            variable_T->getShapeFunctionOrder());
    }
    if (variable_p->getShapeFunctionOrder() != 1)
    {
        OGS_FATAL(
            "The shape function order of pressure must be 1 but its input "
            "value in <process_variable><order> is {:d}. Please correct it.",
            variable_p->getShapeFunctionOrder());
    }

    DBUG("Associate displacement with process variable '{:s}'.",
         variable_u->getName());

    if (variable_u->getNumberOfGlobalComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable_u->getName(),
            variable_u->getNumberOfGlobalComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable '{:s}'.",
         variable_p->getName());
    if (variable_p->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '{:s}' is not a scalar variable but has "
            "{:d} components.",
            variable_p->getName(),
            variable_p->getNumberOfGlobalComponents());
    }

    DBUG("Associate temperature with process variable '{:s}'.",
         variable_T->getName());
    if (variable_T->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "temperature process variable '{:s}' is not a scalar variable but "
            "has {:d} components.",
            variable_T->getName(),
            variable_T->getNumberOfGlobalComponents());
    }

    /// \section parametersthm Process Parameters
    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, materialIDs(mesh), config);

    auto ice_constitutive_relation =
        MaterialLib::Solids::createConstitutiveRelationIce<DisplacementDim>(
            parameters, local_coordinate_system, config);

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_HYDRO_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto const is_linear =
        //! \ogs_file_param{prj__processes__process__linear}
        config.getConfigParameter("linear", false);

    checkThermoOsmosisProperties(media);
    checkLiquidThermalExpansivity(media);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    auto const equation_balance_type_str =
        //! \ogs_file_param{prj__processes__process__THERMO_HYDRO_MECHANICS__equation_balance_type}
        config.getConfigParameter<std::string>("equation_balance_type",
                                               "volume");
    if (equation_balance_type_str == "volume")
    {
        ProcessLib::Common::HydraulicProcess::checkVolumeBalanceEquationSetting(
            media_map);
    }

    // Initial stress conditions

    auto initial_stress = ProcessLib::createInitialStress<DisplacementDim>(
        config, parameters, mesh);

    auto stabilizer = NumLib::createNumericalStabilization(mesh, config);

    bool const is_volume_balance_equation_type =
        equation_balance_type_str == "volume";
    ThermoHydroMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        std::move(ice_constitutive_relation),
        is_volume_balance_equation_type,
        std::move(initial_stress),
        specific_body_force,
        std::move(stabilizer),
        std::move(coupling_scheme),
        thermal_process_id,
        hydraulic_process_id,
        mechanical_process_id};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<ThermoHydroMechanicsProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables), is_linear);
}

template std::unique_ptr<Process> createThermoHydroMechanicsProcess<2>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createThermoHydroMechanicsProcess<3>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);
}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
