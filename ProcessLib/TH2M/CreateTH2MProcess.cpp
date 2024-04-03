/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTH2MProcess.h"

#include <cassert>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Common/HydroMechanics/CreateInitialStress.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/TH2M/ConstitutiveRelations/NoPhaseTransition.h"
#include "ProcessLib/TH2M/ConstitutiveRelations/PhaseTransition.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "TH2MProcess.h"
#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
std::unique_ptr<ConstitutiveRelations::PhaseTransitionModel>
createPhaseTransitionModel(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    // the approach here is that the number of phase components determines the
    // nature of the phase transition: If the gas phase consists of two or more
    // components, evaporation is involved; if the water phase consists of at
    // least two components, gas can be dissolved in water.

    // Fluid phases are always defined in the first medium of the media vector,
    // thus media.begin() points to the right medium.
    const bool phase_transition =
        (media.begin()->second->phase("Gas").numberOfComponents() > 1) &&
        (media.begin()->second->phase("AqueousLiquid").numberOfComponents() >
         1);
    // Only if both fluids consist of more than one component, the model
    // phase_transition is returned.
    if (phase_transition)
    {
        return std::make_unique<ConstitutiveRelations::PhaseTransition>(media);
    }
    return std::make_unique<ConstitutiveRelations::NoPhaseTransition>(media);
}

template <int DisplacementDim>
std::unique_ptr<Process> createTH2MProcess(
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
    config.checkConfigParameter("type", "TH2M");
    DBUG("Create TH2M Process.");
    DBUG(" ");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__TH2M__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

    /// \section processvariablesth2m Process Variables

    //! \ogs_file_param{prj__processes__process__TH2M__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_pGR;
    ProcessVariable* variable_pCap;
    ProcessVariable* variable_T;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        /// Primary process variables as they appear in the global component
        /// vector:
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__gas_pressure}
             "gas_pressure",
             //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__capillary_pressure}
             "capillary_pressure",
             //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__displacement}
             "displacement"});
        variable_pGR = &per_process_variables[0].get();
        variable_pCap = &per_process_variables[1].get();
        variable_T = &per_process_variables[2].get();
        variable_u = &per_process_variables[3].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL("A Staggered version of TH2M is not implemented.");

        using namespace std::string_literals;
        for (auto const& variable_name :
             {"gas_pressure"s, "capillary_pressure"s, "temperature"s,
              "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_pGR = &process_variables[0][0].get();
        variable_pCap = &process_variables[1][0].get();
        variable_T = &process_variables[2][0].get();
        variable_u = &process_variables[3][0].get();
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

    DBUG("Associate gas pressure with process variable '{:s}'.",
         variable_pGR->getName());
    if (variable_pGR->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Gas pressure process variable '{:s}' is not a scalar variable but "
            "has "
            "{:d} components.",
            variable_pGR->getName(),
            variable_pGR->getNumberOfGlobalComponents());
    }

    DBUG("Associate capillary pressure with process variable '{:s}'.",
         variable_pCap->getName());
    if (variable_pCap->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Capillary pressure process variable '{:s}' is not a scalar "
            "variable but has "
            "{:d} components.",
            variable_pCap->getName(),
            variable_pCap->getNumberOfGlobalComponents());
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

    /// \section parametersth2m Process Parameters
    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    // reference temperature
    const auto& reference_temperature = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH2M__reference_temperature}
        "reference_temperature", parameters, 1, &mesh);
    DBUG("Use '{:s}' as reference temperature parameter.",
         reference_temperature.name);

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__TH2M__specific_body_force}
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

    // Initial stress conditions
    auto initial_stress = ProcessLib::createInitialStress<DisplacementDim>(
        config, parameters, mesh);

    auto const mass_lumping =
        //! \ogs_file_param{prj__processes__process__TH2M__mass_lumping}
        config.getConfigParameter<bool>("mass_lumping", false);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    auto phase_transition_model = createPhaseTransitionModel(media);

    const bool use_TaylorHood_elements =
        variable_pCap->getShapeFunctionOrder() !=
                variable_u->getShapeFunctionOrder()
            ? true
            : false;

    TH2MProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        std::move(phase_transition_model),
        reference_temperature,
        std::move(initial_stress),
        specific_body_force,
        mass_lumping,
        use_TaylorHood_elements};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<TH2MProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createTH2MProcess<2>(
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

template std::unique_ptr<Process> createTH2MProcess<3>(
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
}  // namespace TH2M
}  // namespace ProcessLib
