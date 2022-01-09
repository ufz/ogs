/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <optional>

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
namespace PhreeqcIOData
{
struct UserPunch;

std::unique_ptr<UserPunch> createUserPunch(
    std::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
