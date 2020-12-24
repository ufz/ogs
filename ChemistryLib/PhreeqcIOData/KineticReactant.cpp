/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ostream>

#include "KineticReactant.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void KineticReactant::print(std::ostream& os,
                            std::size_t const global_id) const
{
    os << name << "\n";

    if (!chemical_formula.empty())
    {
        os << "-formula " << chemical_formula << "\n";
    }

    os << "-m  " << (*molality)[global_id] << "\n";

    os << "-parms";
    for (auto const& parameter : parameters)
    {
        os << " " << parameter;
    }
    os << "\n";
}
const ItemType KineticReactant::item_type;
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
