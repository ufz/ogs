/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
std::ostream& operator<<(std::ostream& os, UserPunch const& user_punch)
{
    os << "USER_PUNCH" << "\n";
    os << "-headings ";
    auto const& secondary_variables = user_punch.secondary_variables;
    for (auto& secondary_variable : secondary_variables)
        os << secondary_variable.name.c_str() << " ";
    os << "\n";

    os << "-start" << "\n";
    int line_number = 1;
    for (auto const& statement : user_punch.statements)
    {
        os << line_number << " " << statement.c_str() << "\n";
        ++line_number;
    }
    os << "-end" << "\n";

    return os;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
