/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include "KineticReactant.h"

namespace ChemistryLib
{
std::ofstream& operator<<(std::ofstream& out,
                          std::vector<KineticReactant> const& kinetic_reactants)
{
    for (auto const& kinetic_reactant : kinetic_reactants)
    {
        out << kinetic_reactant.name << "\n";

        out << "-m  " << kinetic_reactant.amount << "\n";

        if (!kinetic_reactant.parameters.empty())
        {
            out << "-parms";
            for (auto const& parameter : kinetic_reactant.parameters)
            {
                out << " " << parameter;
            }
            out << "\n";
        }
    }

    return out;
}
}  // namespace ChemistryLib
