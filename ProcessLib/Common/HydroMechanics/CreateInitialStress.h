/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 13, 2024, 10:20 AM
 */

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
