// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ReactionRate.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
ReactionRate::ReactionRate(std::string kinetic_reactant_,
                           std::vector<std::string> const& statements)
    : kinetic_reactant(std::move(kinetic_reactant_))
{
    int line_number = 1;
    for (auto const& statement : statements)
    {
        _commands += std::to_string(line_number) + " " + statement + "; ";
        ++line_number;
    }
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
