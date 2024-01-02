/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
                 std::vector<std::string> const& statements);

    std::string const& commands() const { return _commands; }

public:
    std::string const kinetic_reactant;

private:
    std::string _commands;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
