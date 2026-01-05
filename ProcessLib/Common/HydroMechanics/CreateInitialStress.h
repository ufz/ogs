// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}
namespace ParameterLib
{
struct ParameterBase;
}  // namespace ParameterLib

namespace ProcessLib
{
struct InitialStress;

template <int DisplacementDim>
InitialStress createInitialStress(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh, bool const mandatory_stress_type = false);

extern template InitialStress createInitialStress<2>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh, bool const mandatory_stress_type);
extern template InitialStress createInitialStress<3>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh, bool const mandatory_stress_type);
};  // namespace ProcessLib
