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

#include <iosfwd>
#include <string>
#include <vector>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct ReactionRate
{
    ReactionRate(std::string kinetic_reactant_,
                 std::vector<std::string>&& expression_statements_)
        : kinetic_reactant(std::move(kinetic_reactant_)),
          expression_statements(std::move(expression_statements_))
    {
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    ReactionRate const& reaction_rate);

    std::string const kinetic_reactant;
    std::vector<std::string> const expression_statements;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
