// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
