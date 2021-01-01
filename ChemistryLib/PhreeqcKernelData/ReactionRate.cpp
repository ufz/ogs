/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRate.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
ReactionRate::ReactionRate(std::string kinetic_reactant_,
                           std::vector<std::string> statements)
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
