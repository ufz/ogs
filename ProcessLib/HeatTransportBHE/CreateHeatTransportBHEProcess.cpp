/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHeatTransportBHEProcess.h"

#include "BHE/BHEAbstract.h"
#include "BHE/BHE_1U.h"
#include "BHE/BHE_2U.h"
#include "BHE/BHE_CXA.h"
#include "BHE/BHE_CXC.h"
#include "BHE/CreateBHE1U.h"
#include "BHE/CreateBHE2U.h"
#include "BHE/CreateBHECXA.h"
#include "BHE/CreateBHECXC.h"
#include "BaseLib/Algorithm.h"
#include "HeatTransportBHEProcess.h"
#include "HeatTransportBHEProcessData.h"
#include "MaterialLib/Fluid/Density/CreateFluidDensityModel.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/Fluid/ThermalConductivity/CreateFluidThermalConductivityModel.h"
#include "MaterialLib/Fluid/Viscosity/CreateViscosityModel.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
std::unique_ptr<Process> createHeatTransportBHEProcess(
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
    config.checkConfigParameter("type", "HEAT_TRANSPORT_BHE");

    DBUG("Create HeatTransportBHE Process.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    // reading primary variables for each
    // BHE----------------------------------------------------------
    auto range =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables__process_variable}
        pv_config.getConfigParameterList<std::string>("process_variable");
    std::vector<std::reference_wrapper<ProcessVariable>> per_process_variables;

    for (std::string const& pv_name : range)
    {
        if (pv_name != "temperature_soil" &&
            pv_name.find("temperature_BHE") != 0)
        {
            OGS_FATAL(
                "Found a process variable name '%s'. It should be "
                "'temperature_soil' or 'temperature_BHE_X'");
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

        per_process_variables.emplace_back(
            const_cast<ProcessVariable&>(*variable));
    }
    process_variables.push_back(std::move(per_process_variables));
    // end of reading primary variables for each
    // BHE----------------------------------------------------------

    // solid phase thermal conductivity parameter.
    auto& thermal_conductivity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special
        //! prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_solid}
        "thermal_conductivity_solid", parameters, 1);

    DBUG("Use \'%s\' as solid phase thermal conductivity parameter.",
         thermal_conductivity_solid.name.c_str());

    // solid phase thermal conductivity parameter.
    auto& thermal_conductivity_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_fluid}
        "thermal_conductivity_fluid", parameters, 1);

    DBUG("Use \'%s\' as fluid phase thermal conductivity parameter.",
         thermal_conductivity_fluid.name.c_str());

    // gas phase thermal conductivity parameter.
    auto& thermal_conductivity_gas = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_gas}
        "thermal_conductivity_gas", parameters, 1);

    DBUG("Use \'%s\' as gas phase thermal conductivity parameter.",
         thermal_conductivity_gas.name.c_str());

    // solid phase heat capacity parameter.
    auto& heat_capacity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_solid}
        "heat_capacity_solid", parameters, 1);

    DBUG("Use \'%s\' as solid phase heat capacity parameter.",
         heat_capacity_solid.name.c_str());

    // fluid phase heat capacity parameter.
    auto& heat_capacity_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_fluid}
        "heat_capacity_fluid", parameters, 1);

    DBUG("Use \'%s\' as fluid phase heat capacity parameter.",
         heat_capacity_fluid.name.c_str());

    // gas phase heat capacity parameter.
    auto& heat_capacity_gas = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_gas}
        "heat_capacity_gas", parameters, 1);

    DBUG("Use \'%s\' as gas phase heat capacity parameter.",
         heat_capacity_gas.name.c_str());

    // solid phase density parameter.
    auto& density_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_solid}
        "density_solid", parameters, 1);

    DBUG("Use \'%s\' as solid phase density parameter.",
         density_solid.name.c_str());

    // fluid phase density parameter.
    auto& density_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_fluid}
        "density_fluid", parameters, 1);

    DBUG("Use \'%s\' as fluid phase density parameter.",
         density_fluid.name.c_str());

    // gas phase density parameter.
    auto& density_gas = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_gas}
        "density_gas", parameters, 1);

    DBUG("Use \'%s\' as gas phase density parameter.",
         density_gas.name.c_str());

    // get the refrigerant properties from fluid property class
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_density}
    auto const& rho_conf = fluid_config.getConfigSubtree("refrigerant_density");
    auto bhe_refrigerant_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_viscosity}
    auto const& mu_conf =
        fluid_config.getConfigSubtree("refrigerant_viscosity");
    auto bhe_refrigerant_viscosity =
        MaterialLib::Fluid::createViscosityModel(mu_conf);
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_specific_heat_capacity}
    auto const& cp_conf =
        fluid_config.getConfigSubtree("refrigerant_specific_heat_capacity");
    auto bhe_refrigerant_heat_capacity =
        MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(cp_conf);
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_thermal_conductivity}
    auto const& lambda_conf =
        fluid_config.getConfigSubtree("refrigerant_thermal_conductivity");
    auto bhe_regrigerant_heat_conductivity =
        MaterialLib::Fluid::createFluidThermalConductivityModel(lambda_conf);

    // reading BHE
    // parameters--------------------------------------------------------------
    std::vector<std::unique_ptr<BHE::BHEAbstract>> vec_BHEs;
    // BHE::BHE_Net BHE_network;

    // now read the BHE configurations
    auto const& bhe_configs =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers}
        config.getConfigSubtree("borehole_heat_exchangers");

    for (
        auto const& bhe_conf :
        //! \ogs_file_param
        //! prj__processes__process___HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger}
        bhe_configs.getConfigSubtreeList("borehole_heat_exchanger"))
    {
        // read in the parameters
        using namespace BHE;
        const std::string bhe_type_str =
            bhe_conf.getConfigParameter<std::string>("bhe_type");

        // convert BHE type
        if (bhe_type_str == "BHE_TYPE_1U")
        {
            // initialize the 1U type BHE
            BHE::BHE_1U* m_bhe_1u =
                BHE::CreateBHE1U(bhe_conf,
                                 curves,
                                 bhe_refrigerant_density,
                                 bhe_refrigerant_viscosity,
                                 bhe_refrigerant_heat_capacity,
                                 bhe_regrigerant_heat_conductivity);

            vec_BHEs.emplace_back(std::make_unique<BHE_1U>(*m_bhe_1u));
        }
        else if (bhe_type_str == "BHE_TYPE_2U")
        {
            // initialize the 2U type BHE
            BHE::BHE_2U* m_bhe_2u =
                BHE::CreateBHE2U(bhe_conf,
                                 curves,
                                 bhe_refrigerant_density,
                                 bhe_refrigerant_viscosity,
                                 bhe_refrigerant_heat_capacity,
                                 bhe_regrigerant_heat_conductivity);

            vec_BHEs.emplace_back(std::make_unique<BHE_2U>(*m_bhe_2u));
        }
        else if (bhe_type_str == "BHE_TYPE_CXC")
        {
            // initialize the CXC type BHE
            BHE::BHE_CXC* m_bhe_CXC =
                BHE::CreateBHECXC(bhe_conf,
                                  curves,
                                  bhe_refrigerant_density,
                                  bhe_refrigerant_viscosity,
                                  bhe_refrigerant_heat_capacity,
                                  bhe_regrigerant_heat_conductivity);

            vec_BHEs.emplace_back(std::make_unique<BHE_CXC>(*m_bhe_CXC));
        }
        else if (bhe_type_str == "BHE_TYPE_CXA")
        {
            // initialize the CXA type BHE
            BHE::BHE_CXA* m_bhe_CXA =
                BHE::CreateBHECXA(bhe_conf,
                                  curves,
                                  bhe_refrigerant_density,
                                  bhe_refrigerant_viscosity,
                                  bhe_refrigerant_heat_capacity,
                                  bhe_regrigerant_heat_conductivity);

            vec_BHEs.emplace_back(std::make_unique<BHE_CXA>(*m_bhe_CXA));
        }
    }
    // end of reading BHE
    // parameters-------------------------------------------------------

    HeatTransportBHEProcessData process_data{thermal_conductivity_solid,
                                             thermal_conductivity_fluid,
                                             thermal_conductivity_gas,
                                             heat_capacity_solid,
                                             heat_capacity_fluid,
                                             heat_capacity_gas,
                                             density_solid,
                                             density_fluid,
                                             density_gas,
                                             std::move(vec_BHEs)};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HeatTransportBHE_Temperature"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<HeatTransportBHEProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
