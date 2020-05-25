/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
                 std::vector<std::string> statements);

    std::string const& commands() const { return commands_; }

public:
    std::string const kinetic_reactant;

private:
    std::string commands_;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
