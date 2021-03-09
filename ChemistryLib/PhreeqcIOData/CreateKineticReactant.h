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

#include <optional>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;

template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct KineticReactant;

std::vector<KineticReactant> createKineticReactants(
    std::optional<BaseLib::ConfigTree> const& config, MeshLib::Mesh& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
