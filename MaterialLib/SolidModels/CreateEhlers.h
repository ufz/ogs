/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "CreateNewtonRaphsonSolverParameters.h"
#include "ParameterLib/Utils.h"

#include "Ehlers.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
inline std::unique_ptr<DamagePropertiesParameters> createDamageProperties(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__damage_properties__alpha_d}
    auto& alpha_d =
        ParameterLib::findParameter<double>(config, "alpha_d", parameters, 1);

    DBUG("Use '{:s}' as alpha_d.", alpha_d.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__damage_properties__beta_d}
    auto& beta_d =
        ParameterLib::findParameter<double>(config, "beta_d", parameters, 1);

    DBUG("Use '{:s}' as beta_d.", beta_d.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__damage_properties__h_d}
    auto& h_d =
        ParameterLib::findParameter<double>(config, "h_d", parameters, 1);

    DBUG("Use '{:s}' as h_d.", h_d.name);

    return std::make_unique<DamagePropertiesParameters>(
        DamagePropertiesParameters{alpha_d, beta_d, h_d});
}

template <int DisplacementDim>
std::unique_ptr<SolidEhlers<DisplacementDim>> createEhlers(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "Ehlers");
    DBUG("Create Ehlers material");

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__shear_modulus}
    auto& shear_modulus = ParameterLib::findParameter<double>(
        config, "shear_modulus", parameters, 1);

    DBUG("Use '{:s}' as shear modulus parameter.", shear_modulus.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__bulk_modulus}
    auto& bulk_modulus = ParameterLib::findParameter<double>(
        config, "bulk_modulus", parameters, 1);

    DBUG("Use '{:s}' as bulk modulus parameter.", bulk_modulus.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__kappa}
    auto& kappa =
        ParameterLib::findParameter<double>(config, "kappa", parameters, 1);

    DBUG("Use '{:s}' as kappa.", kappa.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__beta}
    auto& beta =
        ParameterLib::findParameter<double>(config, "beta", parameters, 1);

    DBUG("Use '{:s}' as beta.", beta.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__gamma}
    auto& gamma =
        ParameterLib::findParameter<double>(config, "gamma", parameters, 1);

    DBUG("Use '{:s}' as gamma.", gamma.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__hardening_modulus}
    auto& hardening_modulus = ParameterLib::findParameter<double>(
        config, "hardening_modulus", parameters, 1);

    DBUG("Use '{:s}' as hardening modulus parameter.", hardening_modulus.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__alpha}
    auto& alpha =
        ParameterLib::findParameter<double>(config, "alpha", parameters, 1);

    DBUG("Use '{:s}' as alpha.", alpha.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__delta}
    auto& delta =
        ParameterLib::findParameter<double>(config, "delta", parameters, 1);

    DBUG("Use '{:s}' as delta.", delta.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__eps}
    auto& eps =
        ParameterLib::findParameter<double>(config, "eps", parameters, 1);

    DBUG("Use '{:s}' as eps.", eps.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__m}
    auto& m = ParameterLib::findParameter<double>(config, "m", parameters, 1);

    DBUG("Use '{:s}' as m.", m.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__alphap}
    auto& alphap =
        ParameterLib::findParameter<double>(config, "alphap", parameters, 1);

    DBUG("Use '{:s}' as alphap.", alphap.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__deltap}
    auto& deltap =
        ParameterLib::findParameter<double>(config, "deltap", parameters, 1);

    DBUG("Use '{:s}' as deltap.", deltap.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__epsp}
    auto& epsp =
        ParameterLib::findParameter<double>(config, "epsp", parameters, 1);

    DBUG("Use '{:s}' as epsp.", epsp.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__mp}
    auto& paremeter_mp =
        ParameterLib::findParameter<double>(config, "mp", parameters, 1);

    DBUG("Use '{:s}' as mp.", paremeter_mp.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__betap}
    auto& betap =
        ParameterLib::findParameter<double>(config, "betap", parameters, 1);

    DBUG("Use '{:s}' as betap.", betap.name);

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__gammap}
    auto& gammap =
        ParameterLib::findParameter<double>(config, "gammap", parameters, 1);

    DBUG("Use '{:s}' as gammap.", gammap.name);

    auto tangent_type =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__tangent_type}
        makeTangentType(config.getConfigParameter<std::string>("tangent_type"));

    MaterialPropertiesParameters mp{
        shear_modulus, bulk_modulus, alpha,  beta,
        gamma,         delta,        eps,    m,
        alphap,        betap,        gammap, deltap,
        epsp,          paremeter_mp, kappa,  hardening_modulus};

    // Damage properties.
    std::unique_ptr<DamagePropertiesParameters> ehlers_damage_properties;

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

    return std::make_unique<SolidEhlers<DisplacementDim>>(
        nonlinear_solver_parameters,
        mp,
        std::move(ehlers_damage_properties),
        tangent_type);
}

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
