/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateConstitutiveRelation.h
 *  Created on July 10, 2018, 12:09 PM
 */

#pragma once

#include <map>
#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
struct ParameterBase;
}

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;

template <int DisplacementDim>
std::map<int,
         std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
createConstitutiveRelations(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

extern template std::map<int,
                         std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>>
createConstitutiveRelations<2>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

extern template std::map<int,
                         std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>>
createConstitutiveRelations<3>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);
}
}
