/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CreateTwoPhaseFlowWithPPProcess.h"
#include <cassert>

#include "MaterialLib/TwoPhaseModels/CreateTwoPhaseFlowMaterialProperties.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "TwoPhaseFlowWithPPProcess.h"
#include "TwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
std::unique_ptr<Process> CreateTwoPhaseFlowWithPPProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "TWOPHASE_FLOW_PP");

    DBUG("Create TwoPhaseFlowProcess with PP model.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_PP__process_variables__process_variable}
         "gas_pressure", "capillary_pressure"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"TwoPhaseFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);
    // Specific body force
    Eigen::VectorXd specific_body_force;

    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(b.size() > 0 && b.size() < 4);
    specific_body_force.resize(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
        std::copy_n(b.data(), b.size(), specific_body_force.data());

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__mass_lumping}
    auto mass_lump = config.getConfigParameter<bool>("mass_lumping");

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    auto const& mat_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");

    std::unique_ptr<
        MaterialLib::TwoPhaseFlowWithPP::TwoPhaseFlowWithPPMaterialProperties>
        material = nullptr;

    if (mat_ids)
    {
        INFO("The twophase flow is in heterogeneous porous media.");
        const bool has_material_ids = true;
        material = MaterialLib::TwoPhaseFlowWithPP::
            CreateTwoPhaseFlowMaterialProperties(mat_config, has_material_ids,
                                                 *mat_ids);
        TwoPhaseFlowWithPPProcessData process_data{
            specific_body_force,
            has_gravity,
            mass_lump,
            std::move(material),
            *curves.at("curve_PC_S"),
            *curves.at("curve_S_Krel_wet"),
            *curves.at("curve_S_Krel_nonwet")};
        return std::unique_ptr<Process>{new TwoPhaseFlowWithPPProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller),
            mat_config, curves}};
    }
    else
    {
        INFO("The twophase flow is in homogeneous porous media.");

        MeshLib::Properties dummy_property;

        auto const& dummy_property_vector =
            dummy_property.createNewPropertyVector<int>(
                "MaterialIDs", MeshLib::MeshItemType::Cell, 1);

        // Since dummy_property_vector is only visible in this function,
        // the following constant, has_material_ids, is employed to indicate
        // that material_ids does not exist.
        const bool has_material_ids = false;
        material = MaterialLib::TwoPhaseFlowWithPP::
            CreateTwoPhaseFlowMaterialProperties(mat_config, has_material_ids,
                                                 *dummy_property_vector);
        TwoPhaseFlowWithPPProcessData process_data{
            specific_body_force,
            has_gravity,
            mass_lump,
            std::move(material),
            *curves.at("curve_PC_S"),
            *curves.at("curve_S_Krel_wet"),
            *curves.at("curve_S_Krel_nonwet")};
        return std::unique_ptr<Process>{new TwoPhaseFlowWithPPProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller),
            mat_config, curves}};
    }
}

}  // end of namespace
}  // end of namespace
