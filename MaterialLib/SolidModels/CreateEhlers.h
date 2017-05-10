/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <logog/include/logog.hpp>

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

#include "Ehlers.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
inline NumLib::NewtonRaphsonSolverParameters
createNewtonRaphsonSolverParameters(BaseLib::ConfigTree const& config)
{
    DBUG("Create local nonlinear solver parameters.");
    auto const maximum_iterations =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__nonlinear_solver__maximum_iterations}
        config.getConfigParameter<int>("maximum_iterations");

    DBUG("\tmaximum_iterations: %d.", maximum_iterations);

    auto const error_tolerance =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__nonlinear_solver__error_tolerance}
        config.getConfigParameter<double>("error_tolerance");

    DBUG("\terror_tolerance: %g.", error_tolerance);

    return {maximum_iterations, error_tolerance};
}

inline std::unique_ptr<EhlersDamageProperties> createDamageProperties(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__damage_properties__alpha_d}
    auto& alpha_d =
        ProcessLib::findParameter<double>(config, "alpha_d", parameters, 1);

    DBUG("Use \'%s\' as alpha_d.", alpha_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__damage_properties__beta_d}
    auto& beta_d =
        ProcessLib::findParameter<double>(config, "beta_d", parameters, 1);

    DBUG("Use \'%s\' as beta_d.", beta_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__damage_properties__h_d}
    auto& h_d = ProcessLib::findParameter<double>(config, "h_d", parameters, 1);

    DBUG("Use \'%s\' as h_d.", h_d.name.c_str());

    return std::make_unique<EhlersDamageProperties>(
        EhlersDamageProperties{alpha_d, beta_d, h_d});
}

template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createEhlers(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "Ehlers");
    DBUG("Create Ehlers material");

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__shear_modulus}
    auto& shear_modulus = ProcessLib::findParameter<double>(
        config, "shear_modulus", parameters, 1);

    DBUG("Use \'%s\' as shear modulus parameter.", shear_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__bulk_modulus}
    auto& bulk_modulus = ProcessLib::findParameter<double>(
        config, "bulk_modulus", parameters, 1);

    DBUG("Use \'%s\' as bulk modulus parameter.", bulk_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__kappa}
    auto& kappa =
        ProcessLib::findParameter<double>(config, "kappa", parameters, 1);

    DBUG("Use \'%s\' as kappa.", kappa.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__beta}
    auto& beta =
        ProcessLib::findParameter<double>(config, "beta", parameters, 1);

    DBUG("Use \'%s\' as beta.", beta.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__gamma}
    auto& gamma =
        ProcessLib::findParameter<double>(config, "gamma", parameters, 1);

    DBUG("Use \'%s\' as gamma.", gamma.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__hardening_modulus}
    auto& hardening_modulus = ProcessLib::findParameter<double>(
        config, "hardening_modulus", parameters, 1);

    DBUG("Use \'%s\' as hardening modulus parameter.",
         hardening_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__alpha}
    auto& alpha =
        ProcessLib::findParameter<double>(config, "alpha", parameters, 1);

    DBUG("Use \'%s\' as alpha.", alpha.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__delta}
    auto& delta =
        ProcessLib::findParameter<double>(config, "delta", parameters, 1);

    DBUG("Use \'%s\' as delta.", delta.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__eps}
    auto& eps = ProcessLib::findParameter<double>(config, "eps", parameters, 1);

    DBUG("Use \'%s\' as eps.", eps.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__m}
    auto& m = ProcessLib::findParameter<double>(config, "m", parameters, 1);

    DBUG("Use \'%s\' as m.", m.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__alphap}
    auto& alphap =
        ProcessLib::findParameter<double>(config, "alphap", parameters, 1);

    DBUG("Use \'%s\' as alphap.", alphap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__deltap}
    auto& deltap =
        ProcessLib::findParameter<double>(config, "deltap", parameters, 1);

    DBUG("Use \'%s\' as deltap.", deltap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__epsp}
    auto& epsp =
        ProcessLib::findParameter<double>(config, "epsp", parameters, 1);

    DBUG("Use \'%s\' as epsp.", epsp.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__mp}
    auto& paremeter_mp =
        ProcessLib::findParameter<double>(config, "mp", parameters, 1);

    DBUG("Use \'%s\' as mp.", paremeter_mp.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__betap}
    auto& betap =
        ProcessLib::findParameter<double>(config, "betap", parameters, 1);

    DBUG("Use \'%s\' as betap.", betap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__gammap}
    auto& gammap =
        ProcessLib::findParameter<double>(config, "gammap", parameters, 1);

    DBUG("Use \'%s\' as gammap.", gammap.name.c_str());

    MaterialProperties mp{
        shear_modulus, bulk_modulus, alpha,  beta,
        gamma,         delta,        eps,    m,
        alphap,        betap,        gammap, deltap,
        epsp,          paremeter_mp, kappa,  hardening_modulus};

    // Damage properties.
    std::unique_ptr<EhlersDamageProperties> ehlers_damage_properties;

    auto const& ehlers_damage_config =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__damage_properties}
        config.getConfigSubtreeOptional("damage_properties");
    if (ehlers_damage_config)
    {
        ehlers_damage_properties =
            createDamageProperties(parameters, *ehlers_damage_config);
    }

    auto const& nonlinear_solver_config =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__nonlinear_solver}
        config.getConfigSubtree("nonlinear_solver");
    auto const nonlinear_solver_parameters =
        createNewtonRaphsonSolverParameters(nonlinear_solver_config);

    return std::unique_ptr<MechanicsBase<DisplacementDim>>{
        new SolidEhlers<DisplacementDim>{nonlinear_solver_parameters, mp,
                                         std::move(ehlers_damage_properties)}};
}

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
