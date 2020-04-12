/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "CreateNewtonRaphsonSolverParameters.h"
#include "Lubby2.h"
#include "ParameterLib/Utils.h"

namespace MaterialLib
{
namespace Solids
{
namespace Lubby2
{
template <int DisplacementDim>
std::unique_ptr<Lubby2<DisplacementDim>> createLubby2(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "Lubby2");
    DBUG("Create Lubby2 material");

    // Kelvin shear modulus.
    auto& kelvin_shear_modulus = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__kelvin_shear_modulus}
        config, "kelvin_shear_modulus", parameters, 1);

    DBUG("Use '{:s}' as kelvin shear modulus parameter.",
         kelvin_shear_modulus.name);

    // Kelvin viscosity.
    auto& kelvin_viscosity = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__kelvin_viscosity}
        config, "kelvin_viscosity", parameters, 1);

    DBUG("Use '{:s}' as kelvin viscosity parameter.", kelvin_viscosity.name);

    // Maxwell shear modulus.
    auto& maxwell_shear_modulus = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__maxwell_shear_modulus}
        config, "maxwell_shear_modulus", parameters, 1);

    DBUG("Use '{:s}' as maxwell shear modulus parameter.",
         maxwell_shear_modulus.name);

    // Maxwell bulk modulus.
    auto& maxwell_bulk_modulus = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__maxwell_bulk_modulus}
        config, "maxwell_bulk_modulus", parameters, 1);

    DBUG("Use '{:s}' as maxwell bulk modulus parameter.",
         maxwell_bulk_modulus.name);

    // Maxwell viscosity.
    auto& maxwell_viscosity = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__maxwell_viscosity}
        config, "maxwell_viscosity", parameters, 1);

    DBUG("Use '{:s}' as maxwell viscosity parameter.", maxwell_viscosity.name);

    // Dependency parameter for mK.
    auto& dependency_parameter_mK = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mk}
        config, "dependency_parameter_mk", parameters, 1);

    DBUG("Use '{:s}' as dependency parameter mK.",
         dependency_parameter_mK.name);

    // Dependency parameter for mvK.
    auto& dependency_parameter_mvK = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mvk}
        config, "dependency_parameter_mvk", parameters, 1);

    DBUG("Use '{:s}' as dependency parameter mvK.",
         dependency_parameter_mvK.name);

    // Dependency parameter for mvM.
    auto& dependency_parameter_mvM = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mvm}
        config, "dependency_parameter_mvm", parameters, 1);

    DBUG("Use '{:s}' as dependency parameter mvM.",
         dependency_parameter_mvM.name);

    auto const dependency_parameter_tref_name = config.getConfigParameterOptional<
        std::string>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_tref}
        "dependency_parameter_tref");

    auto const dependency_parameter_mgt_name = config.getConfigParameterOptional<
        std::string>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mgt}
        "dependency_parameter_mgt");

    auto const dependency_parameter_mkt_name = config.getConfigParameterOptional<
        std::string>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mkt}
        "dependency_parameter_mkt");

    auto const dependency_parameter_q_name = config.getConfigParameterOptional<
        std::string>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_q}
        "dependency_parameter_q");

    Lubby2MaterialProperties* mp;

    // Check all the thermal related parameter settings
    if (dependency_parameter_tref_name || dependency_parameter_mgt_name ||
        dependency_parameter_mkt_name || dependency_parameter_q_name)
    {
        WARN("Lubby2 model lost essential thermal parameter(s), the thermal related properties are not activated.");
    }

    // Initialize the thermal related parameters 
    if (dependency_parameter_tref_name && dependency_parameter_mgt_name &&
        dependency_parameter_mkt_name && dependency_parameter_q_name)
    {
        // Dependency parameter for Tref
        auto& dependency_parameter_Tref = ParameterLib::findParameter<double>(
            *dependency_parameter_tref_name, parameters, 1);

        DBUG("Use '%s' as dependency parameter Tref.",
             dependency_parameter_Tref.name.c_str());

        // Dependency parameter for mGT
        auto& dependency_parameter_mGT = ParameterLib::findParameter<double>(
            *dependency_parameter_mgt_name, parameters, 1);

        DBUG("Use '%s' as dependency parameter mGT.",
             dependency_parameter_mGT.name.c_str());

        // Dependency parameter for mKT
        auto& dependency_parameter_mKT = ParameterLib::findParameter<double>(
            *dependency_parameter_mkt_name, parameters, 1);

        DBUG("Use '%s' as dependency parameter mKT.",
             dependency_parameter_mKT.name.c_str());

        // Dependency parameter for Q
        auto& dependency_parameter_Q = ParameterLib::findParameter<double>(
            *dependency_parameter_q_name, parameters, 1);

        DBUG("Use '%s' as dependency parameter Q.",
             dependency_parameter_Q.name.c_str());

        Lubby2ThermalMaterialProperties mp_thermal{
            dependency_parameter_Tref, dependency_parameter_mGT,
            dependency_parameter_mKT, dependency_parameter_Q};

        mp = new Lubby2MaterialProperties(kelvin_shear_modulus,
                                          maxwell_shear_modulus,
                                          maxwell_bulk_modulus,
                                          kelvin_viscosity,
                                          maxwell_viscosity,
                                          dependency_parameter_mK,
                                          dependency_parameter_mvK,
                                          dependency_parameter_mvM,
                                          mp_thermal);
    }
    else
    {
        mp = new Lubby2MaterialProperties(kelvin_shear_modulus,
                                          maxwell_shear_modulus,
                                          maxwell_bulk_modulus,
                                          kelvin_viscosity,
                                          maxwell_viscosity,
                                          dependency_parameter_mK,
                                          dependency_parameter_mvK,
                                          dependency_parameter_mvM);
    }

    auto const& nonlinear_solver_config =
        //! \ogs_file_param{material__solid__constitutive_relation__Lubby2__nonlinear_solver}
        config.getConfigSubtree("nonlinear_solver");
    auto const nonlinear_solver_parameters =
        createNewtonRaphsonSolverParameters(nonlinear_solver_config);

    return std::unique_ptr<Lubby2<DisplacementDim>>{
        new Lubby2<DisplacementDim>{nonlinear_solver_parameters, *mp}};
}

}  // namespace Lubby2
}  // namespace Solids
}  // namespace MaterialLib
