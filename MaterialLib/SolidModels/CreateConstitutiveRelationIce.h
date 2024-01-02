/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

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
}  // namespace ParameterLib

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;

template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createConstitutiveRelationIce(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>
createConstitutiveRelationIce(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>
createConstitutiveRelationIce(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace Solids
}  // namespace MaterialLib
