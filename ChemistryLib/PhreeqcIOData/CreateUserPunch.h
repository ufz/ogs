// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
