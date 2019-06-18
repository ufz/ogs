/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ostream>

#include "EquilibriumPhase.h"

namespace ChemistryLib
{
std::ostream& operator<<(std::ostream& os,
                         EquilibriumPhase const& equilibrium_phase)
{
    os << equilibrium_phase.name;

    os << " " << equilibrium_phase.saturation_index;

    os << " " << equilibrium_phase.amount << "\n";

    return os;
}
}  // namespace ChemistryLib
