// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
