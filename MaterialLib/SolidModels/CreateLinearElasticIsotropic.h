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

#include "ParameterLib/Utils.h"

#include "LinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<LinearElasticIsotropic<DisplacementDim>>
createLinearElasticIsotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config, const bool skip_type_checking)
{
    if (!skip_type_checking)
    {
        //! \ogs_file_param{material__solid__constitutive_relation__type}
        config.checkConfigParameter("type", "LinearElasticIsotropic");
        DBUG("Create LinearElasticIsotropic material");
    }

    // Youngs modulus
    auto& youngs_modulus = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticIsotropic__youngs_modulus}
        config, "youngs_modulus", parameters, 1);

    DBUG("Use '{:s}' as youngs_modulus parameter.", youngs_modulus.name);

    // Poissons ratio
    auto& poissons_ratio = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticIsotropic__poissons_ratio}
        config, "poissons_ratio", parameters, 1);

    DBUG("Use '{:s}' as poissons_ratio parameter.", poissons_ratio.name);

    typename LinearElasticIsotropic<DisplacementDim>::MaterialProperties mp{
        youngs_modulus, poissons_ratio};

    return std::make_unique<LinearElasticIsotropic<DisplacementDim>>(mp);
}

}  // namespace Solids
}  // namespace MaterialLib
