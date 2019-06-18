/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ostream>

#include "KineticReactant.h"

namespace ChemistryLib
{
std::ostream& operator<<(std::ostream& os,
                         KineticReactant const& kinetic_reactant)
{
    os << kinetic_reactant.name << "\n";

    if (!kinetic_reactant.chemical_formula.empty())
    {
        os << "-formula " << kinetic_reactant.chemical_formula << "\n";
    }

    os << "-m  " << kinetic_reactant.amount << "\n";

    if (!kinetic_reactant.parameters.empty())
    {
        os << "-parms";
        for (auto const& parameter : kinetic_reactant.parameters)
        {
            os << " " << parameter;
        }
        os << "\n";
    }

    return os;
}
}  // namespace ChemistryLib
