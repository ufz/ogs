/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHydroMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
std::unique_ptr<Process> createHydroMechanicsProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HYDRO_MECHANICS");
    DBUG("Create HydroMechanicsProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__process_variables__pressure}
            "pressure",
            //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__process_variables__displacement}
            "displacement"});
        variable_p = &per_process_variables[0].get();
        variable_u = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_p = &process_variables[0][0].get();
        variable_u = &process_variables[1][0].get();
    }

    DBUG("Associate displacement with process variable '%s'.",
         variable_u->getName().c_str());

    if (variable_u->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            variable_u->getName().c_str(),
            variable_u->getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable '%s'.",
         variable_p->getName().c_str());
    if (variable_p->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_p->getName().c_str(),
            variable_p->getNumberOfComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    // Intrinsic permeability
    auto& intrinsic_permeability = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1, &mesh);

    DBUG("Use '%s' as intrinsic conductivity parameter.",
         intrinsic_permeability.name.c_str());

    // Fluid density
    auto& fluid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__fluid_density}
        "fluid_density", parameters, 1, &mesh);
    DBUG("Use '%s' as fluid density parameter.", fluid_density.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    std::array const requiredGasProperties = {MaterialPropertyLib::viscosity};
    std::array const requiredSolidProperties = {
        MaterialPropertyLib::porosity, MaterialPropertyLib::biot_coefficient,
        MaterialPropertyLib::density};
    for (auto const& m : media)
    {
        m.second->phase("Gas").checkRequiredProperties(requiredGasProperties);
        m.second->phase("Solid").checkRequiredProperties(
            requiredSolidProperties);
    }

    // Initial stress conditions
    auto const initial_stress = ParameterLib::findOptionalTagParameter<double>(
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__initial_stress}
        config, "initial_stress", parameters,
        // Symmetric tensor size, 4 or 6, not a Kelvin vector.
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value,
        &mesh);

    // Reference temperature
    double const reference_temperature =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__reference_temperature}
        config.getConfigParameter<double>(
            "reference_temperature", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'reference_temperature' as reference temperature.");

    //Specific gas constant
    double const specific_gas_constant =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__specific_gas_constant}
        config.getConfigParameter<double>(
            "specific_gas_constant", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'specific_gas_constant' as specific gas constant.");

    // Fluid compressibility
    double const fluid_compressibility =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__fluid_compressibility}
        config.getConfigParameter<double>(
            "fluid_compressibility", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'fluid_compressibility' as fluid compressibility parameter.");

    auto const fluid_type =
        FluidType::strToFluidType(
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__fluid_type}
            config.getConfigParameter<std::string>("fluid_type"));
    DBUG("Use 'fluid_type' as fluid type parameter.");

    if (!FluidType::checkRequiredParams(fluid_type, fluid_compressibility,
                                        reference_temperature,
                                        specific_gas_constant))
    {
        OGS_FATAL(FluidType::getErrorMsg(fluid_type));
    }

    HydroMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),     std::move(media_map),
        std::move(solid_constitutive_relations),
        initial_stress,        intrinsic_permeability,
        fluid_density,         specific_body_force,
        fluid_compressibility, reference_temperature,
        specific_gas_constant, fluid_type};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HydroMechanicsProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createHydroMechanicsProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createHydroMechanicsProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace HydroMechanics
}  // namespace ProcessLib
