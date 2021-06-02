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
