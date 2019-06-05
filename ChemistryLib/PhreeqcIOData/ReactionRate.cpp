/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include "ReactionRate.h"

namespace ChemistryLib
{
std::ofstream& operator<<(std::ofstream& out,
                          std::vector<ReactionRate> const& reaction_rates)
{
    for (auto const& reaction_rate : reaction_rates)
    {
        out << reaction_rate.kinetic_reactant << "\n";
        out << "-start" << "\n";
        int line_number = 1;
        for (auto const& expression_statement :
             reaction_rate.expression_statements)
        {
            out << line_number << " " << expression_statement << "\n";
            ++line_number;
        }
        out << "-end" << "\n";
    }

    return out;
}
}  // namespace ChemistryLib
