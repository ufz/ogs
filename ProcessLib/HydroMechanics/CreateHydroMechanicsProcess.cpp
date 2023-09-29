/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHydroMechanicsProcess.h"

#include <cassert>
#include <string>

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcessData.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MeshLib/Utils/Is2DMeshOnRotatedVerticalPlane.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace HydroMechanics
{

CouplingScheme parseCouplingScheme(
    std::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return Monolithic{};
    }

    auto const coupling_scheme_type =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__type}
        config->getConfigParameter<std::string>("type");

    if (coupling_scheme_type == "monolithic")
    {
        return Monolithic{};
    }

    // Default value is 0.5, which is recommended by [Mikelic & Wheeler].
    double const fixed_stress_stabilization_parameter =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__fixed_stress_stabilization_parameter}
        config->getConfigParameter<double>(
            "fixed_stress_stabilization_parameter", 0.5);

    DBUG("Using value {:g} for coupling parameter of staggered scheme.",
         fixed_stress_stabilization_parameter);

    {  // Check parameter value. Optimum is not a-priori known, but within
       // certain interval [Storvik & Nordbotten]
        double const csp_min = 1.0 / 6.0;
        double const csp_max = 1.0;
        if (fixed_stress_stabilization_parameter < csp_min ||
            fixed_stress_stabilization_parameter > csp_max)
        {
            WARN(
                "Value of coupling scheme parameter = {:g} is out of "
                "reasonable range ({:g}, {:g}).",
                fixed_stress_stabilization_parameter, csp_min, csp_max);
        }
    }

    bool const fixed_stress_over_time_step =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__fixed_stress_over_time_step}
        config->getConfigParameter<std::string>("fixed_stress_over_time_step",
                                                "false") == "true";

    return Staggered{fixed_stress_stabilization_parameter,
                     fixed_stress_over_time_step};
}

template <int DisplacementDim>
std::unique_ptr<Process> createHydroMechanicsProcess(
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
    config.checkConfigParameter("type", "HYDRO_MECHANICS");
    DBUG("Create HydroMechanicsProcess.");

    if (DisplacementDim == 2)
    {
        if (mesh.isAxiallySymmetric() &&
            MeshLib::is2DMeshOnRotatedVerticalPlane(mesh))
        {
            OGS_FATAL(
                "Mesh {:s} is on a plane rotated around the vertical axis. The "
                "axisymmetric problem can not use such mesh.",
                mesh.getName());
        }
    }

    auto const coupling_scheme = parseCouplingScheme(
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme}
        config.getConfigSubtreeOptional("coupling_scheme"));

    /// \section processvariableshm Process Variables

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p = nullptr;
    ProcessVariable* variable_u = nullptr;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    int const hydraulic_process_id = 0;
    int mechanics_related_process_id = 0;

    if (std::holds_alternative<Monolithic>(coupling_scheme))
    {
        /// Primary process variables as they appear in the global component
        /// vector:
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

    if (std::holds_alternative<Staggered>(coupling_scheme))
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        mechanics_related_process_id = 1;
        variable_p = &process_variables[hydraulic_process_id][0].get();
        variable_u = &process_variables[mechanics_related_process_id][0].get();
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

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    /// \section parametershm Process Parameters
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
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__mass_lumping}
    auto mass_lumping = config.getConfigParameter<bool>("mass_lumping", false);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    std::array const requiredMediumProperties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::permeability, MaterialPropertyLib::porosity,
        MaterialPropertyLib::biot_coefficient};
    std::array const requiredFluidProperties = {MaterialPropertyLib::viscosity,
                                                MaterialPropertyLib::density};
    std::array const requiredSolidProperties = {MaterialPropertyLib::density};

    for (auto const& element_id : mesh.getElements() | MeshLib::views::ids)
    {
        media_map.checkElementHasMedium(element_id);
        auto const& medium = *media_map.getMedium(element_id);
        checkRequiredProperties(medium, requiredMediumProperties);
        checkRequiredProperties(fluidPhase(medium), requiredFluidProperties);
        checkRequiredProperties(medium.phase("Solid"), requiredSolidProperties);
    }
    DBUG("Media properties verified.");

    // Initial stress conditions
    auto const initial_stress = ParameterLib::findOptionalTagParameter<double>(
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__initial_stress}
        config, "initial_stress", parameters,
        // Symmetric tensor size, 4 or 6, not a Kelvin vector.
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim),
        &mesh);

    const bool use_taylor_hood_elements = variable_p->getShapeFunctionOrder() !=
                                          variable_u->getShapeFunctionOrder();

    HydroMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        initial_stress,
        specific_body_force,
        mass_lumping,
        coupling_scheme,
        hydraulic_process_id,
        mechanics_related_process_id,
        use_taylor_hood_elements};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HydroMechanicsProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::holds_alternative<Monolithic>(coupling_scheme));
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
}  // namespace ProcessLib
