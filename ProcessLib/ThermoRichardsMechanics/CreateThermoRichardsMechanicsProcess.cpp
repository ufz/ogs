/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoRichardsMechanicsProcess.h"

#include <cassert>

#include "ConstitutiveStress_StrainTemperature/CreateConstitutiveSetting.h"
#include "ConstitutiveStress_StrainTemperature/Traits.h"

#ifdef OGS_USE_MFRONT
#include "ConstitutiveStressSaturation_StrainPressureTemperature/CreateConstitutiveSetting.h"
#include "ConstitutiveStressSaturation_StrainPressureTemperature/Traits.h"
#endif

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Common/HydroMechanics/CreateInitialStress.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "ThermoRichardsMechanicsProcess.h"
#include "ThermoRichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::porosity, MaterialPropertyLib::biot_coefficient,
        MaterialPropertyLib::bishops_effective_stress,
        MaterialPropertyLib::relative_permeability,
        MaterialPropertyLib::saturation};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::viscosity, MaterialPropertyLib::density};
    std::array const required_solid_properties = {MaterialPropertyLib::density};

    // Thermal properties are not checked because they can be phase property or
    // medium property (will be enabled later).
    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, required_medium_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
        checkRequiredProperties(m.second->phase("Solid"),
                                required_solid_properties);
    }
}

void checkProcessVariableComponents(ProcessVariable const& variable,
                                    const int dim)
{
    DBUG("Associate displacement with process variable '{:s}'.",
         variable.getName());

    if (variable.getNumberOfGlobalComponents() != dim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable.getName(),
            variable.getNumberOfGlobalComponents(),
            dim);
    }
}

template <int DisplacementDim, typename ConstitutiveTraits,
          typename CreateConstitutiveSetting>
std::unique_ptr<Process> createThermoRichardsMechanicsProcessStage2(
    std::string const& name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    bool const mandatory_stress0_type)
{
    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

    /// \section processvariablestrm Process Variables

    //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_T;
    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        /// Primary process variables as they appear in the global component
        /// vector:
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__THERMO_RICHARDS_MECHANICS__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__THERMO_RICHARDS_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__THERMO_RICHARDS_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_T = &per_process_variables[0].get();
        variable_p = &per_process_variables[1].get();
        variable_u = &per_process_variables[2].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL(
            "So far, only the monolithic scheme is implemented for "
            "THERMO_RICHARDS_MECHANICS");
    }

    checkProcessVariableComponents(*variable_T, 1);
    checkProcessVariableComponents(*variable_p, 1);
    checkProcessVariableComponents(*variable_u, DisplacementDim);

    /// \section parameterstrm Process Parameters

    auto solid_constitutive_relations =
        CreateConstitutiveSetting::createSolidConstitutiveRelations(
            parameters, local_coordinate_system, config);

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__specific_body_force}
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

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    DBUG(
        "Check the media properties of ThermoRichardsMechanics process "
        "...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    auto initial_stress = ProcessLib::createInitialStress<DisplacementDim>(
        config, parameters, mesh, mandatory_stress0_type);

    bool mass_lumping = false;
    if (auto const mass_lumping_ptr =
            //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__mass_lumping}
        config.getConfigParameterOptional<bool>("mass_lumping"))
    {
        DBUG("Using mass lumping for the Richards flow equation.");
        mass_lumping = *mass_lumping_ptr;
    }

    bool const apply_body_force_for_deformation =
        //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__apply_body_force_for_deformation}
        config.getConfigParameter<bool>("apply_body_force_for_deformation",
                                        true);

    bool const initialize_porosity_from_medium_property =
        //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__initialize_porosity_from_medium_property}
        config.getConfigParameter("initialize_porosity_from_medium_property",
                                  true);

    const bool use_TaylorHood_elements =
        variable_p->getShapeFunctionOrder() !=
                variable_u->getShapeFunctionOrder()
            ? true
            : false;

    ThermoRichardsMechanicsProcessData<DisplacementDim, ConstitutiveTraits>
        process_data{materialIDs(mesh),
                     std::move(media_map),
                     std::move(solid_constitutive_relations),
                     std::move(initial_stress),
                     specific_body_force,
                     mass_lumping,
                     use_TaylorHood_elements,
                     apply_body_force_for_deformation,
                     InitializePorosityFromMediumProperty{
                         initialize_porosity_from_medium_property}};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<
        ThermoRichardsMechanicsProcess<DisplacementDim, ConstitutiveTraits>>(
        name, mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

template <int DisplacementDim>
std::unique_ptr<Process> createThermoRichardsMechanicsProcess(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_RICHARDS_MECHANICS");
    DBUG("Create ThermoRichardsMechanicsProcess.");

    auto const subtype =
        //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__subtype}
        config.getConfigParameter<std::string>("subtype",
                                               "Stress_StrainTemperature");
    INFO("TRM process subtype is '{}'", subtype);

    if (subtype == "Stress_StrainTemperature")
    {
        bool const mandatory_stress0_type = false;
        return createThermoRichardsMechanicsProcessStage2<
            DisplacementDim,
            ConstitutiveStress_StrainTemperature::ConstitutiveTraits<
                DisplacementDim>,
            ConstitutiveStress_StrainTemperature::CreateConstitutiveSetting<
                DisplacementDim>>(name, mesh, std::move(jacobian_assembler),
                                  variables, parameters,
                                  local_coordinate_system, integration_order,
                                  config, media, mandatory_stress0_type);
    }

    if (subtype == "StressSaturation_StrainPressureTemperature")
    {
#ifdef OGS_USE_MFRONT
        bool const mandatory_stress0_type = true;
        return createThermoRichardsMechanicsProcessStage2<
            DisplacementDim,
            ConstitutiveStressSaturation_StrainPressureTemperature::
                ConstitutiveTraits<DisplacementDim>,
            ConstitutiveStressSaturation_StrainPressureTemperature::
                CreateConstitutiveSetting<DisplacementDim>>(
            name, mesh, std::move(jacobian_assembler), variables, parameters,
            local_coordinate_system, integration_order, config, media,
            mandatory_stress0_type);
#else
        OGS_FATAL(
            "TRM process subtype 'StressSaturation_StrainPressureTemperature' "
            "is not supported, because OGS has not been built with MFront.");
#endif
    }

    OGS_FATAL("Unknown TRM process subtype '{}'.", subtype);
}

template std::unique_ptr<Process> createThermoRichardsMechanicsProcess<2>(
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

template std::unique_ptr<Process> createThermoRichardsMechanicsProcess<3>(
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

}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
