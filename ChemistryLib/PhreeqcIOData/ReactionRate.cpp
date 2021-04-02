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

#include <ostream>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::ostream& operator<<(std::ostream& os, ReactionRate const& reaction_rate)
{
    os << reaction_rate.kinetic_reactant << "\n";
    os << "-start"
       << "\n";
    int line_number = 1;
    for (auto const& expression_statement : reaction_rate.expression_statements)
    {
        os << line_number << " " << expression_statement << "\n";
        ++line_number;
    }
    os << "-end"
       << "\n";

    return os;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
