/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHydroMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/FractureModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/CreateMohrCoulomb.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <unsigned GlobalDim>
std::unique_ptr<Process> createHydroMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HYDRO_MECHANICS_WITH_LIE");
    DBUG("Create HydroMechanicsProcess with LIE.");
    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        (staggered_scheme && (*staggered_scheme == "staggered")) ? false : true;

    // Process variables
    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__process_variables}
    auto const pv_conf = config.getConfigSubtree("process_variables");
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
            pv_name.find("displacement_jump") != 0)
        {
            OGS_FATAL(
                "Found a process variable name '%s'. It should be "
                "'displacement' or 'displacement_jumpN' or 'pressure'");
        }
        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [&pv_name](ProcessVariable const& v) {
                                         return v.getName() == pv_name;
                                     });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '%s' in the provided "
                "variables "
                "list for config tag <%s>.",
                pv_name.c_str(), "process_variable");
        }
        DBUG("Found process variable \'%s\' for config tag <%s>.",
             variable->getName().c_str(), "process_variable");

        if (pv_name.find("displacement") != std::string::npos &&
            variable->getNumberOfComponents() != GlobalDim)
        {
            OGS_FATAL(
                "Number of components of the process variable '%s' is "
                "different "
                "from the displacement dimension: got %d, expected %d",
                variable->getName().c_str(),
                variable->getNumberOfComponents(),
                GlobalDim);
        }

        if (!use_monolithic_scheme)
        {
            if (pv_name == "pressure")
                p_process_variables.emplace_back(
                    const_cast<ProcessVariable&>(*variable));
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
        OGS_FATAL("Currently only one displacement jump is supported");

    if (!use_monolithic_scheme)
    {
        process_variables.push_back(std::move(p_process_variables));
        process_variables.push_back(std::move(u_process_variables));
    }
    else
        process_variables.push_back(std::move(p_u_process_variables));


    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>> material =
        nullptr;
    if (type == "LinearElasticIsotropic")
    {
        material = MaterialLib::Solids::createLinearElasticIsotropic<GlobalDim>(
            parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Intrinsic permeability
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1);

    DBUG("Use \'%s\' as intrinsic permeabiltiy parameter.",
         intrinsic_permeability.name.c_str());

    // Storage coefficient
    auto& specific_storage = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__specific_storage}
        "specific_storage", parameters, 1);

    DBUG("Use \'%s\' as specific storage parameter.",
         specific_storage.name.c_str());

    // Fluid viscosity
    auto& fluid_viscosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fluid_viscosity}
        "fluid_viscosity", parameters, 1);
    DBUG("Use \'%s\' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Fluid density
    auto& fluid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fluid_density}
        "fluid_density", parameters, 1);
    DBUG("Use \'%s\' as fluid density parameter.", fluid_density.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__biot_coefficient}
        "biot_coefficient", parameters, 1);
    DBUG("Use \'%s\' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Porosity
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__porosity}
        "porosity", parameters, 1);
    DBUG("Use \'%s\' as porosity parameter.", porosity.name.c_str());

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__solid_density}
        "solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.", solid_density.name.c_str());

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
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
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
        else if (frac_type == "MohrCoulomb")
        {
            fracture_model =
                MaterialLib::Fracture::createMohrCoulomb<GlobalDim>(
                    parameters, fracture_model_config);
        }
        else
        {
            OGS_FATAL(
                "Cannot construct fracture constitutive relation of given type "
                "\'%s\'.",
                frac_type.c_str());
        }
    }

    // Fracture properties
    std::unique_ptr<FractureProperty> frac_prop = nullptr;
    auto opt_fracture_properties_config =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties}
        config.getConfigSubtreeOptional("fracture_properties");
    if (opt_fracture_properties_config)
    {
        auto& fracture_properties_config = *opt_fracture_properties_config;

        frac_prop = std::make_unique<ProcessLib::LIE::FractureProperty>();
        frac_prop->mat_id =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__material_id}
            fracture_properties_config.getConfigParameter<int>("material_id");
        frac_prop->aperture0 = &ProcessLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__initial_aperture}
            fracture_properties_config, "initial_aperture", parameters, 1);
        frac_prop->specific_storage = &ProcessLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__specific_storage}
            fracture_properties_config, "specific_storage", parameters, 1);
        frac_prop->biot_coefficient = &ProcessLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__biot_coefficient}
            fracture_properties_config, "biot_coefficient", parameters, 1);
    }

    // initial effective stress in matrix
    auto& initial_effective_stress = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__initial_effective_stress}
        "initial_effective_stress", parameters,
        ProcessLib::KelvinVectorDimensions<GlobalDim>::value);
    DBUG("Use \'%s\' as initial effective stress parameter.",
         initial_effective_stress.name.c_str());

    // initial effective stress in fracture
    auto& initial_fracture_effective_stress = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__initial_fracture_effective_stress}
        "initial_fracture_effective_stress", parameters, GlobalDim);
    DBUG("Use \'%s\' as initial fracture effective stress parameter.",
         initial_fracture_effective_stress.name.c_str());

    // deactivation of matrix elements in flow
    auto opt_deactivate_matrix_in_flow =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__deactivate_matrix_in_flow}
        config.getConfigParameterOptional<bool>("deactivate_matrix_in_flow");
    bool const deactivate_matrix_in_flow =
        opt_deactivate_matrix_in_flow && *opt_deactivate_matrix_in_flow;
    ;
    if (deactivate_matrix_in_flow)
        INFO("Deactivate matrix elements in flow calculation.");

    HydroMechanicsProcessData<GlobalDim> process_data{
        std::move(material),
        intrinsic_permeability,
        specific_storage,
        fluid_viscosity,
        fluid_density,
        biot_coefficient,
        porosity,
        solid_density,
        specific_body_force,
        std::move(fracture_model),
        std::move(frac_prop),
        initial_effective_stress,
        initial_fracture_effective_stress,
        deactivate_matrix_in_flow};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HydroMechanics_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::make_unique<HydroMechanicsProcess<GlobalDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createHydroMechanicsProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);
template std::unique_ptr<Process> createHydroMechanicsProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
