/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "CreateNewtonRaphsonSolverParameters.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

#include "Lubby2.h"

namespace MaterialLib
{
namespace Solids
{
namespace Lubby2
{
template <int DisplacementDim>
std::unique_ptr<Lubby2<DisplacementDim>> createLubby2(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "Lubby2");
    DBUG("Create Lubby2 material");

    // Kelvin shear modulus.
    auto& kelvin_shear_modulus = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__kelvin_shear_modulus}
        config, "kelvin_shear_modulus", parameters, 1);

    DBUG("Use '%s' as kelvin shear modulus parameter.",
         kelvin_shear_modulus.name.c_str());

    // Kelvin viscosity.
    auto& kelvin_viscosity = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__kelvin_viscosity}
        config, "kelvin_viscosity", parameters, 1);

    DBUG("Use '%s' as kelvin viscosity parameter.",
         kelvin_viscosity.name.c_str());

    // Maxwell shear modulus.
    auto& maxwell_shear_modulus = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__maxwell_shear_modulus}
        config, "maxwell_shear_modulus", parameters, 1);

    DBUG("Use '%s' as maxwell shear modulus parameter.",
         maxwell_shear_modulus.name.c_str());

    // Maxwell bulk modulus.
    auto& maxwell_bulk_modulus = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__maxwell_bulk_modulus}
        config, "maxwell_bulk_modulus", parameters, 1);

    DBUG("Use '%s' as maxwell bulk modulus parameter.",
         maxwell_bulk_modulus.name.c_str());

    // Maxwell viscosity.
    auto& maxwell_viscosity = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__maxwell_viscosity}
        config, "maxwell_viscosity", parameters, 1);

    DBUG("Use '%s' as maxwell viscosity parameter.",
         maxwell_viscosity.name.c_str());

    // Dependency parameter for mK.
    auto& dependency_parameter_mK = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mk}
        config, "dependency_parameter_mk", parameters, 1);

    DBUG("Use '%s' as dependency parameter mK.",
         dependency_parameter_mK.name.c_str());

    // Dependency parameter for mvK.
    auto& dependency_parameter_mvK = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mvk}
        config, "dependency_parameter_mvk", parameters, 1);

    DBUG("Use '%s' as dependency parameter mvK.",
         dependency_parameter_mvK.name.c_str());

    // Dependency parameter for mvM.
    auto& dependency_parameter_mvM = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__Lubby2__dependency_parameter_mvm}
        config, "dependency_parameter_mvm", parameters, 1);

    DBUG("Use '%s' as dependency parameter mvM.",
         dependency_parameter_mvM.name.c_str());

    Lubby2MaterialProperties mp{
        kelvin_shear_modulus,     maxwell_shear_modulus,
        maxwell_bulk_modulus,     kelvin_viscosity,
        maxwell_viscosity,        dependency_parameter_mK,
        dependency_parameter_mvK, dependency_parameter_mvM};

    auto const nonlinear_solver_parameters =
        createNewtonRaphsonSolverParameters(config);

    return std::unique_ptr<Lubby2<DisplacementDim>>{
        new Lubby2<DisplacementDim>{nonlinear_solver_parameters, mp}};
}

}  // namespace Lubby2
}  // namespace Solids
}  // namespace MaterialLib
