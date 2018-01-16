/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSmallDeformationProcess.h"

#include <cassert>

#include "MaterialLib/FractureModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/CreateMohrCoulomb.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"

#include "ProcessLib/Output/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

#include "SmallDeformationProcess.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <int DisplacementDim>
std::unique_ptr<Process> createSmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "SMALL_DEFORMATION_WITH_LIE");
    DBUG("Create SmallDeformationProcess with LIE.");

    // Process variables
    //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__process_variables}
    auto const pv_conf = config.getConfigSubtree("process_variables");
    auto range =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__process_variables__process_variable}
        pv_conf.getConfigParameterList<std::string>("process_variable");
    std::vector<std::reference_wrapper<ProcessVariable>> per_process_variables;

    for (std::string const& pv_name : range)
    {
        if (pv_name != "displacement" && pv_name.find("displacement_jump") != 0)
        {
            OGS_FATAL(
                "Found a process variable name '%s'. It should be "
                "'displacement' or 'displacement_jumpN'");
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

        per_process_variables.emplace_back(const_cast<ProcessVariable&>(*variable));
    }
    auto const n_fractures = per_process_variables.size() - 1;
    if (n_fractures < 1)
    {
        OGS_FATAL("No displacement jump variables are specified");
    }

    DBUG("Associate displacement with process variable \'%s\'.",
         per_process_variables.back().get().getName().c_str());

    if (per_process_variables.back().get().getNumberOfComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            per_process_variables.back().get().getName().c_str(),
            per_process_variables.back().get().getNumberOfComponents(),
            DisplacementDim);
    }
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Fracture constitutive relation.
    // read type;
    auto const fracture_model_config =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__fracture_model}
        config.getConfigSubtree("fracture_model");

    auto const frac_type =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__fracture_model__type}
        fracture_model_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<DisplacementDim>>
        fracture_model = nullptr;
    if (frac_type == "LinearElasticIsotropic")
    {
        fracture_model = MaterialLib::Fracture::createLinearElasticIsotropic<
            DisplacementDim>(parameters, fracture_model_config);
    }
    else if (frac_type == "MohrCoulomb")
    {
        fracture_model =
            MaterialLib::Fracture::createMohrCoulomb<DisplacementDim>(
                parameters, fracture_model_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct fracture constitutive relation of given type "
            "\'%s\'.",
            frac_type.c_str());
    }

    // Fracture properties
    std::vector<std::unique_ptr<FractureProperty>> vec_fracture_property;
    for (
        auto fracture_properties_config :
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__fracture_properties}
        config.getConfigSubtreeList("fracture_properties"))
    {
        auto& para_b0 = ProcessLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__fracture_properties__initial_aperture}
            fracture_properties_config, "initial_aperture", parameters, 1);
        auto frac_prop(new FractureProperty());
        frac_prop->fracture_id = vec_fracture_property.size();
        frac_prop->mat_id =
            //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__fracture_properties__material_id}
            fracture_properties_config.getConfigParameter<int>("material_id");
        frac_prop->aperture0 = &para_b0;
        vec_fracture_property.emplace_back(frac_prop);
    }

    if (n_fractures != vec_fracture_property.size())
    {
        OGS_FATAL(
            "The number of displacement jumps and the number of "
            "<fracture_properties> "
            "are not consistent");
    }

    SmallDeformationProcessData<DisplacementDim> process_data(
        std::move(material), std::move(fracture_model),
        std::move(vec_fracture_property));

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"SmallDeformation_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::make_unique<SmallDeformationProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}

template std::unique_ptr<Process> createSmallDeformationProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createSmallDeformationProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
