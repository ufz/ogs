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

#include <memory>
#include <optional>

#include "KineticReactant.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
std::unique_ptr<Kinetics> createKineticReactants(
    std::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh);
}
}  // namespace ChemistryLib
