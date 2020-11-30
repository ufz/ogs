/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UserPunch.h"
#include "BaseLib/ConfigTreeUtil.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void UserPunch::initialize(std::size_t const num_chemical_systems)
{
    for (auto& secondary_variable : secondary_variables)
    {
        secondary_variable.value->resize(num_chemical_systems);
    }
}

std::ostream& operator<<(std::ostream& os, UserPunch const& user_punch)
{
    os << "USER_PUNCH" << "\n";
    os << "-headings ";
    auto const& secondary_variables = user_punch.secondary_variables;
    for (auto& secondary_variable : secondary_variables)
    {
        os << secondary_variable.name << " ";
    }
    os << "\n";

    os << "-start" << "\n";
    int line_number = 1;
    for (auto const& statement : user_punch.statements)
    {
        os << line_number << " " << statement << "\n";
        ++line_number;
    }
    os << "-end" << "\n";

    return os;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
