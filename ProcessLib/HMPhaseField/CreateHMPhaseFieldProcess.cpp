/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHMPhaseFieldProcess.h"

#include <cassert>

#include "HMPhaseFieldProcess.h"
#include "HMPhaseFieldProcessData.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace HMPhaseField
{
void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const requiredMediumProperties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::porosity, MaterialPropertyLib::permeability};
    std::array const requiredFluidProperties = {MaterialPropertyLib::viscosity,
                                                MaterialPropertyLib::density};
    std::array const requiredSolidProperties = {
        MaterialPropertyLib::biot_coefficient, MaterialPropertyLib::density};

    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, requiredMediumProperties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                requiredFluidProperties);
        checkRequiredProperties(m.second->phase("Solid"),
                                requiredSolidProperties);
    }
}

template <int DisplacementDim>
std::unique_ptr<Process> createHMPhaseFieldProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HM_PHASE_FIELD");
    DBUG("Create HMPhaseFieldProcess.");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__coupling_scheme}
        config.getConfigParameter<std::string>("coupling_scheme", "staggered");
    const bool use_monolithic_scheme = (coupling_scheme != "staggered");

    //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_ph;
    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        OGS_FATAL("Monolithic implementation is not available.");
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (
            /// Primary process variables as they appear in the global component
            /// vector:
            auto const& variable_name :
            {//! \ogs_file_param_special{prj__processes__process__HM_PHASE_FIELD__process_variables__phasefield}
             "phasefield"s,
             //! \ogs_file_param_special{prj__processes__process__HM_PHASE_FIELD__process_variables__pressure}
             "pressure"s,
             //! \ogs_file_param_special{prj__processes__process__HM_PHASE_FIELD__process_variables__displacement}
             "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_ph = &process_variables[0][0].get();
        variable_p = &process_variables[1][0].get();
        variable_u = &process_variables[2][0].get();
    }

    DBUG("Associate phase field with process variable '{}'.",
         variable_ph->getName());
    if (variable_ph->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "HMPhasefield process variable '{}' is not a scalar variable but "
            "has {:d} components.",
            variable_ph->getName(),
            variable_ph->getNumberOfGlobalComponents());
    }

    DBUG("Associate pressure with process variable '{}'.",
         variable_p->getName());
    if (variable_p->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "HMPhasefield process variable '{}' is not a scalar variable but "
            "has {:d} components.",
            variable_p->getName(),
            variable_p->getNumberOfGlobalComponents());
    }

    DBUG("Associate displacement with process variable '{}'.",
         variable_u->getName());

    if (variable_u->getNumberOfGlobalComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "The component number of the process variable '{}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable_u->getName(),
            variable_u->getNumberOfGlobalComponents(),
            DisplacementDim);
    }

    /// \section parametershmpf Process Parameters
    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, materialIDs(mesh), config);

    auto const phasefield_parameters_config =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__phasefield_parameters}
        config.getConfigSubtree("phasefield_parameters");

    // Residual stiffness
    auto const& residual_stiffness = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__HM_PHASE_FIELD__phasefield_parameters__residual_stiffness}
        "residual_stiffness", parameters, 1);
    DBUG("Use '{}' as residual stiffness.", residual_stiffness.name);

    // Crack resistance
    auto const& crack_resistance = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__HM_PHASE_FIELD__phasefield_parameters__crack_resistance}
        "crack_resistance", parameters, 1);
    DBUG("Use '{}' as crack resistance.", crack_resistance.name);

    // Crack length scale
    auto const& crack_length_scale = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__HM_PHASE_FIELD__phasefield_parameters__crack_length_scale}
        "crack_length_scale", parameters, 1);
    DBUG("Use '{}' as crack length scale.", crack_length_scale.name);

    // Characteristic_length
    auto const characteristic_length =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__characteristic_length}
        config.getConfigParameter<double>("characteristic_length", 1.0);

    // Fluid compression modulus. Default value is 0.0.
    auto const fluid_compressibility =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__fluid_compressibility}
        config.getConfigParameter<double>("fluid_compressibility", 0.0);

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__specific_body_force}
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

    // Specific fracture direction
    Eigen::Vector<double, DisplacementDim> specific_fracture_direction;
    {
        std::vector<double> const fracture_normal =
            //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__specific_fracture_direction}
            config.getConfigParameter<std::vector<double>>(
                "specific_fracture_direction");
        if (fracture_normal.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific fracture direction vector does not "
                "match the displacement dimension. Vector size is {:d}, "
                "displacement dimension is {:d}",
                fracture_normal.size(), DisplacementDim);
        }

        std::copy_n(fracture_normal.data(), fracture_normal.size(),
                    specific_fracture_direction.data());
    }

    auto const irreversible_threshold =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__irreversible_threshold}
        config.getConfigParameter<double>("irreversible_threshold", 0.05);

    auto const phasefield_model_string =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__phasefield_model}
        config.getConfigParameter<std::string>("phasefield_model");
    auto const phasefield_model =
        MaterialLib::Solids::Phasefield::convertStringToPhaseFieldModel<
            DisplacementDim>(phasefield_model_string);

    auto const softening_curve_string =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__softening_curve}
        config.getConfigParameterOptional<std::string>("softening_curve");
    auto const softening_curve =
        MaterialLib::Solids::Phasefield::convertStringToSofteningCurve<
            DisplacementDim>(softening_curve_string);

    auto const energy_split_model_string =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__energy_split_model}
        config.getConfigParameter<std::string>("energy_split_model");
    auto const energy_split_model =
        MaterialLib::Solids::Phasefield::convertStringToEnergySplitModel<
            DisplacementDim>(energy_split_model_string);

    auto degradation_derivative = creatDegradationDerivative<DisplacementDim>(
        phasefield_model, characteristic_length, softening_curve);

    auto const diffused_range_parameter =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__diffused_range_parameter}
        config.getConfigParameter<double>("diffused_range_parameter", 2.0);

    auto const fracture_threshold =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__fracture_threshold}
        config.getConfigParameter<double>("fracture_threshold", 1.0);

    auto const fracture_permeability_parameter =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__fracture_permeability_parameter}
        config.getConfigParameter<double>("fracture_permeability_parameter",
                                          0.0);

    // Default value is 0.5, which is recommended by [Mikelic & Wheeler].
    double const fixed_stress_stabilization_parameter =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__fixed_stress_stabilization_parameter}
        config.getConfigParameter<double>(
            "fixed_stress_stabilization_parameter", 0.5);

    DBUG("Using value {:g} for coupling parameter of staggered scheme.",
         fixed_stress_stabilization_parameter);

    auto spatial_stabilization_parameter_optional =
        //! \ogs_file_param{prj__processes__process__HM_PHASE_FIELD__spatial_stabilization_parameter}
        config.getConfigParameterOptional<double>(
            "spatial_stabilization_parameter");
    double const spatial_stabilization_parameter =
        spatial_stabilization_parameter_optional
            ? spatial_stabilization_parameter_optional.value()
            : 0.0;
    DBUG("Using value {} for spatial stablization coupling parameter.",
         spatial_stabilization_parameter);

    // Initial width
    auto const& width_init = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HM_PHASE_FIELD__width_init}
        "width_init", parameters, 0.0);
    DBUG("Use '{}' as initial width.", width_init.name);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    DBUG("Check the media properties of HMPhaseField process ...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    HMPhaseFieldProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        residual_stiffness,
        crack_resistance,
        crack_length_scale,
        specific_body_force,
        specific_fracture_direction,
        irreversible_threshold,
        phasefield_model,
        energy_split_model,
        softening_curve,
        characteristic_length,
        std::move(degradation_derivative),
        diffused_range_parameter,
        fluid_compressibility,
        fracture_threshold,
        fracture_permeability_parameter,
        fixed_stress_stabilization_parameter,
        spatial_stabilization_parameter,
        width_init};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HMPhaseFieldProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createHMPhaseFieldProcess<2>(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createHMPhaseFieldProcess<3>(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace HMPhaseField
}  // namespace ProcessLib
