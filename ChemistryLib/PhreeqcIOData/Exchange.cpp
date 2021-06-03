/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
std::ostream& operator<<(std::ostream& os, ExchangeSite const& exchange_site)
{
    os << exchange_site.ion_exchanging_species << " " << exchange_site.molality << "\n";

    return os;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
