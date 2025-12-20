// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class ReactionRate
{
public:
    ReactionRate(std::string kinetic_reactant_,
                 std::vector<std::string> const& statements);

    std::string const& commands() const { return _commands; }

public:
    std::string const kinetic_reactant;

private:
    std::string _commands;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
