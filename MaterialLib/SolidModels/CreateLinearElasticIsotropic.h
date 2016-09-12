/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_SOLIDMODELS_CREATELINEARELASTICISOTROPIC_H_
#define MATERIALLIB_SOLIDMODELS_CREATELINEARELASTICISOTROPIC_H_

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter
#include "LinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Solids
{

template <int DisplacementDim>
std::unique_ptr<LinearElasticIsotropic<DisplacementDim>>
createLinearElasticIsotropic(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "LinearElasticIsotropic");
    DBUG("Create LinearElasticIsotropic material");

    // Youngs modulus
    auto& youngs_modulus = ProcessLib::findParameter<double>(
        config, "youngs_modulus", parameters, 1);

    DBUG("Use '%s' as youngs_modulus parameter.", youngs_modulus.name.c_str());

    // Poissons ratio
    auto& poissons_ratio = ProcessLib::findParameter<double>(
        config, "poissons_ratio", parameters, 1);

    DBUG("Use '%s' as poissons_ratio parameter.", poissons_ratio.name.c_str());

    typename LinearElasticIsotropic<DisplacementDim>::MaterialProperties mp{
        youngs_modulus, poissons_ratio};

    return std::unique_ptr<LinearElasticIsotropic<DisplacementDim>>{
        new LinearElasticIsotropic<DisplacementDim>{mp}};
}

}  // namespace Solids
}  // namespace MaterialLib

#endif  // MATERIALLIB_SOLIDMODELS_CREATELINEARELASTICISOTROPIC_H_
