/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHydroMechanicsProcess.h"

#include <cassert>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/transform.hpp>

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcessData.h"
#include "MaterialLib/FractureModels/CreateCohesiveZoneModeI.h"
#include "MaterialLib/FractureModels/CreateCoulomb.h"
#include "MaterialLib/FractureModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/Permeability/CreatePermeabilityModel.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/PropertyType.h"
#include "MaterialLib/MPL/Utils/CheckMPLPhasesForSinglePhaseFlow.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "ParameterLib/SpatialPosition.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <int GlobalDim>
std::unique_ptr<Process> createHydroMechanicsProcess(
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
    config.checkConfigParameter("type", "HYDRO_MECHANICS_WITH_LIE");
    DBUG("Create HydroMechanicsProcess with LIE.");
    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

    /// \section processvariablesliehm Process Variables
    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__process_variables}
    auto const pv_conf = config.getConfigSubtree("process_variables");
    /// Primary process variables as they appear in the global component vector:
    auto range =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__process_variables__process_variable}
        pv_conf.getConfigParameterList<std::string>("process_variable");
    std::vector<std::reference_wrapper<ProcessVariable>> p_u_process_variables;
    std::vector<std::reference_wrapper<ProcessVariable>> p_process_variables;
    std::vector<std::reference_wrapper<ProcessVariable>> u_process_variables;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    for (std::string const& pv_name : range)
    {
        if (pv_name != "pressure" && pv_name != "displacement" &&
            !pv_name.starts_with("displacement_jump"))
        {
            OGS_FATAL(
                "Found a process variable name '{}'. It should be "
                "'displacement' or 'displacement_jumpN' or 'pressure'",
                pv_name);
        }
        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [&pv_name](ProcessVariable const& v)
                                     { return v.getName() == pv_name; });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '{:s}' in the provided "
                "variables list for config tag <{:s}>.",
                pv_name, "process_variable");
        }
        DBUG("Found process variable '{:s}' for config tag <{:s}>.",
             variable->getName(), "process_variable");

        if (pv_name.find("displacement") != std::string::npos &&
            variable->getNumberOfGlobalComponents() != GlobalDim)
        {
            OGS_FATAL(
                "Number of components of the process variable '{:s}' is "
                "different from the displacement dimension: got {:d}, expected "
                "{:d}",
                variable->getName(),
                variable->getNumberOfGlobalComponents(),
                GlobalDim);
        }

        if (!use_monolithic_scheme)
        {
            if (pv_name == "pressure")
            {
                p_process_variables.emplace_back(
                    const_cast<ProcessVariable&>(*variable));
            }
            else
            {
                u_process_variables.emplace_back(
                    const_cast<ProcessVariable&>(*variable));
            }
        }
        else
        {
            p_u_process_variables.emplace_back(
                const_cast<ProcessVariable&>(*variable));
        }
    }

    if (p_u_process_variables.size() > 3 || u_process_variables.size() > 2)
    {
        OGS_FATAL("Currently only one displacement jump is supported");
    }

    if (!use_monolithic_scheme)
    {
        process_variables.push_back(std::move(p_process_variables));
        process_variables.push_back(std::move(u_process_variables));
    }
    else
    {
        process_variables.push_back(std::move(p_u_process_variables));
    }

    /// \section parametersliehm Process Parameters
    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<GlobalDim>(
            parameters, local_coordinate_system, materialIDs(mesh), config);

    // Specific body force
    Eigen::Matrix<double, GlobalDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != GlobalDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), GlobalDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    // Fracture constitutive relation.
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>
        fracture_model = nullptr;
    auto const opt_fracture_model_config =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_model}
        config.getConfigSubtreeOptional("fracture_model");
    if (opt_fracture_model_config)
    {
        auto& fracture_model_config = *opt_fracture_model_config;

        auto const frac_type =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_model__type}
            fracture_model_config.peekConfigParameter<std::string>("type");

        if (frac_type == "LinearElasticIsotropic")
        {
            fracture_model =
                MaterialLib::Fracture::createLinearElasticIsotropic<GlobalDim>(
                    parameters, fracture_model_config);
        }
        else if (frac_type == "Coulomb")
        {
            fracture_model = MaterialLib::Fracture::createCoulomb<GlobalDim>(
                parameters, fracture_model_config);
        }
        else if (frac_type == "CohesiveZoneModeI")
        {
            fracture_model = MaterialLib::Fracture::CohesiveZoneModeI::
                createCohesiveZoneModeI<GlobalDim>(parameters,
                                                   fracture_model_config);
        }
        else
        {
            OGS_FATAL(
                "Cannot construct fracture constitutive relation of given type "
                "'{:s}'.",
                frac_type);
        }
    }

    // Fracture properties
    std::vector<FractureProperty> fracture_properties;
    for (
        auto fracture_properties_config :
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties}
        config.getConfigSubtreeList("fracture_properties"))
    {
        fracture_properties.emplace_back(
            fracture_properties.size(),
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__material_id}
            fracture_properties_config.getConfigParameter<int>("material_id"),
            ParameterLib::findParameter<double>(
                //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__initial_aperture}
                fracture_properties_config, "initial_aperture", parameters, 1,
                &mesh));
    }

    std::size_t const n_var_du =
        ranges::accumulate(
            process_variables | ranges::views::transform(ranges::size),
            std::size_t{0}) -
        2 /* for pressure and displacement */;
    if (n_var_du != fracture_properties.size())
    {
        OGS_FATAL(
            "The number of displacement jumps {} and the number of "
            "<fracture_properties> {} are not consistent.",
            n_var_du,
            fracture_properties.size());
    }

    // initial effective stress in matrix
    auto& initial_effective_stress = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__initial_effective_stress}
        "initial_effective_stress", parameters,
        MathLib::KelvinVector::kelvin_vector_dimensions(GlobalDim), &mesh);
    DBUG("Use '{:s}' as initial effective stress parameter.",
         initial_effective_stress.name);

    // initial effective stress in fracture
    auto& initial_fracture_effective_stress = ParameterLib::findParameter<
        double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__initial_fracture_effective_stress}
        "initial_fracture_effective_stress", parameters, GlobalDim, &mesh);
    DBUG("Use '{:s}' as initial fracture effective stress parameter.",
         initial_fracture_effective_stress.name);

    // deactivation of matrix elements in flow
    auto opt_deactivate_matrix_in_flow =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__deactivate_matrix_in_flow}
        config.getConfigParameterOptional<bool>("deactivate_matrix_in_flow");
    bool const deactivate_matrix_in_flow =
        opt_deactivate_matrix_in_flow && *opt_deactivate_matrix_in_flow;

    if (deactivate_matrix_in_flow)
        INFO("Deactivate matrix elements in flow calculation.");

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__use_b_bar}
    auto const use_b_bar = config.getConfigParameter<bool>("use_b_bar", false);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    std::array const requiredMediumProperties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::permeability,
        MaterialPropertyLib::biot_coefficient};
    std::array const requiredFluidProperties = {MaterialPropertyLib::viscosity,
                                                MaterialPropertyLib::density};
    std::array const requiredSolidProperties = {MaterialPropertyLib::density};

    MaterialPropertyLib::checkMPLPhasesForSinglePhaseFlow(mesh, media_map);

    for (auto const& medium : media_map.media())
    {
        checkRequiredProperties(*medium, requiredMediumProperties);
        checkRequiredProperties(fluidPhase(*medium), requiredFluidProperties);
        checkRequiredProperties(medium->phase("Solid"),
                                requiredSolidProperties);
    }

    // Check whether fracture permeability is given as a scalar value.
    for (auto const& element_id : mesh.getElements() | MeshLib::views::ids)
    {
        media_map.checkElementHasMedium(element_id);
        auto const& medium = *media_map.getMedium(element_id);

        // For fracture element
        if (mesh.getElement(element_id)->getDimension() != GlobalDim)
        {
            ParameterLib::SpatialPosition x_position;
            MaterialPropertyLib::VariableArray variables;
            auto const permeability =
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(variables, x_position, 0.0 /*t*/, 0.0 /*dt*/);
            if (!std::holds_alternative<double>(permeability))
            {
                OGS_FATAL(
                    "The permeability model for the fracture must be "
                    "isotropic, and it must return a scalar value.");
            }
        }
    }

    HydroMechanicsProcessData<GlobalDim> process_data{
        materialIDs(mesh),         std::move(solid_constitutive_relations),
        std::move(media_map),      specific_body_force,
        std::move(fracture_model), std::move(fracture_properties),
        initial_effective_stress,  initial_fracture_effective_stress,
        deactivate_matrix_in_flow, use_b_bar};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HydroMechanicsProcess<GlobalDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createHydroMechanicsProcess<2>(
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

template std::unique_ptr<Process> createHydroMechanicsProcess<3>(
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

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
