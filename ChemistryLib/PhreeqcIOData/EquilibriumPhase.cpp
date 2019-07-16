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
void EquilibriumPhase::print(std::ostream& os,
                             std::size_t const chemical_system_id) const
{
    os << name << " " << saturation_index << " "
       << (*amount)[chemical_system_id] << "\n";
}
}  // namespace ChemistryLib
