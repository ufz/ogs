/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSmallDeformationProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/CreateMohrCoulomb.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

#include "SmallDeformationProcess.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{
template <int DisplacementDim>
class SmallDeformationProcess;

extern template class SmallDeformationProcess<2>;

template <int DisplacementDim>
std::unique_ptr<Process>
createSmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "SMALL_DEFORMATION_WITH_LIE");
    DBUG("Create SmallDeformationProcess with LIE.");

    // Process variables
    auto const pv_conf = config.getConfigSubtree("process_variables");
    auto range = pv_conf.getConfigParameterList<std::string>("process_variable");
    std::vector<std::reference_wrapper<ProcessVariable>> process_variables;
    for (std::string const& pv_name : range)
    {
        if (pv_name != "displacement"
            && pv_name.find("displacement_jump")==std::string::npos)
            OGS_FATAL("Found a process variable name '%s'. It should be 'displacement' or 'displacement_jumpN'");
        auto variable = std::find_if(
            variables.cbegin(), variables.cend(),
            [&pv_name](ProcessVariable const& v) { return v.getName() == pv_name; });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '%s' in the provided variables "
                "list for config tag <%s>.",
                pv_name.c_str(), "process_variable");
        }
        DBUG("Found process variable \'%s\' for config tag <%s>.",
             variable->getName().c_str(), "process_variable");

        process_variables.emplace_back(const_cast<ProcessVariable&>(*variable));
    }
    auto const n_fractures = process_variables.size()-1;
    if (n_fractures < 1)
        OGS_FATAL("No displacement jump variables are specified");

    DBUG("Associate displacement with process variable \'%s\'.",
         process_variables.back().get().getName().c_str());

    if (process_variables.back().get().getNumberOfComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            process_variables.back().get().getName().c_str(),
            process_variables.back().get().getNumberOfComponents(),
            DisplacementDim);
    }

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{process__SMALL_DEFORMATION_WITH_LIE__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>> material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        material = MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
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
    auto const fracture_constitutive_relation_config =
        //! \ogs_file_param{process__SMALL_DEFORMATION_WITH_LIE__constitutive_relation}
        config.getConfigSubtree("fracture_constitutive_relation");

    auto const frac_type =
        fracture_constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<DisplacementDim>> fracture_model = nullptr;
    if (frac_type == "LinearElasticIsotropic")
    {
        fracture_model = MaterialLib::Fracture::createLinearElasticIsotropic<DisplacementDim>(
            parameters, fracture_constitutive_relation_config);
    }
    else if (frac_type == "MohrCoulomb")
    {
        fracture_model = MaterialLib::Fracture::createMohrCoulomb<DisplacementDim>(
            parameters, fracture_constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct fracture constitutive relation of given type \'%s\'.",
            frac_type.c_str());
    }

    // Fracture properties
    //! \ogs_file_param{process__SMALL_DEFORMATION_WITH_LIE__fracture_properties}
    auto fracture_properties_config = config.getConfigSubtree("fracture_properties");
    auto &para_b0 = ProcessLib::findParameter<double>(fracture_properties_config, "initial_aperture", parameters, 1);
    std::unique_ptr<FractureProperty> frac_prop(new FractureProperty());
    frac_prop->mat_id = fracture_properties_config.getConfigParameter<int>("material_id");
    frac_prop->aperture0 = &para_b0;


    SmallDeformationProcessData<DisplacementDim> process_data(
        std::move(material), std::move(fracture_model), std::move(frac_prop));

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"SmallDeformation_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<SmallDeformationProcess<DisplacementDim>>{
        new SmallDeformationProcess<DisplacementDim>{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller)}};
}

template
std::unique_ptr<Process>
createSmallDeformationProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

