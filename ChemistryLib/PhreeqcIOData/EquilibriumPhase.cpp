/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include "EquilibriumPhase.h"

namespace ChemistryLib
{
std::ofstream& operator<<(
    std::ofstream& out, std::vector<EquilibriumPhase> const& equilibrium_phases)
{
    for (auto const& equilibrium_phase : equilibrium_phases)
    {
        out << equilibrium_phase.name;

        out << " " << equilibrium_phase.saturation_index;

        out << " " << equilibrium_phase.amount << "\n";
    }

    return out;
}
}  // namespace ChemistryLib
