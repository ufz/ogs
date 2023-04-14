/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>

namespace MeshLib
{
class Mesh;
}

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
class NumericalStabilization;
}

namespace NumLib
{
std::unique_ptr<NumericalStabilization> createNumericalStabilization(
    MeshLib::Mesh const& mesh, BaseLib::ConfigTree const& config);
}  // namespace NumLib
