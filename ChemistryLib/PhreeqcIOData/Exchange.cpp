/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Exchange.h"

#include <ostream>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void ExchangeSite::print(std::ostream& os,
                         std::size_t const chemical_system_id) const
{
    os << name << " " << (*molality)[chemical_system_id] << "\n";
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
