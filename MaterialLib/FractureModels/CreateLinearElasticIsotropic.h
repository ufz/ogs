/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "FractureModelBase.h"

namespace  MaterialLib
{
namespace Fracture
{

template <int DisplacementDim>
std::unique_ptr<FractureModelBase<DisplacementDim>>
createLinearElasticIsotropic(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace Fracture
}  // namespace MaterialLib
