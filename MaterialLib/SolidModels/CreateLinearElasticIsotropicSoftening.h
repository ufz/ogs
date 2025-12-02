/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LinearElasticIsotropicSoftening.h"
#include "ParameterLib/Utils.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<LinearElasticIsotropicSoftening<DisplacementDim>>
createLinearElasticIsotropicSoftening(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config, const bool skip_type_checking)
{
    if (!skip_type_checking)
    {
        //! \ogs_file_param{material__solid__constitutive_relation__type}
        config.checkConfigParameter("type", "LinearElasticIsotropicSoftening");
        DBUG("Create LinearElasticIsotropic material");
    }

    // Youngs modulus
    auto& youngs_modulus = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticIsotropicSoftening__youngs_modulus}
        config, "youngs_modulus", parameters, 1);

    DBUG("Use '{:s}' as youngs_modulus parameter.", youngs_modulus.name);

    // Poissons ratio
    auto& poissons_ratio = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticIsotropicSoftening__poissons_ratio}
        config, "poissons_ratio", parameters, 1);

    DBUG("Use '{:s}' as poissons_ratio parameter.", poissons_ratio.name);

    // Material strength factor
    auto& strength = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticIsotropicSoftening__strength}
        config, "strength", parameters, 1);

    DBUG("Use '{:s}' as strength parameter.", strength.name);

    typename LinearElasticIsotropicSoftening<
        DisplacementDim>::MaterialProperties mp{youngs_modulus, poissons_ratio};

    return std::make_unique<LinearElasticIsotropicSoftening<DisplacementDim>>(
        mp, strength);
}

}  // namespace Solids
}  // namespace MaterialLib
