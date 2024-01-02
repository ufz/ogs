/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on July 10, 2018, 12:09 PM
 */

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ParameterLib
{
struct ParameterBase;
struct CoordinateSystem;
}

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;

template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

template <int DisplacementDim>
std::map<int,
         std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
createConstitutiveRelations(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::map<int,
                         std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>>
createConstitutiveRelations<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::map<int,
                         std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>>
createConstitutiveRelations<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace Solids
}  // namespace MaterialLib
