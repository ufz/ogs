// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
#include <vector>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Output;
struct ChemicalSystem;
struct UserPunch;

std::unique_ptr<Output> createOutput(
    ChemicalSystem const& chemical_system,
    std::unique_ptr<UserPunch> const& user_punch,
    bool const use_high_precision,
    std::string const& project_file_name);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
