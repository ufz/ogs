// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EquilibriumReactant.h"

#include <ostream>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void EquilibriumReactant::print(std::ostream& os,
                                std::size_t const chemical_system_id) const
{
    os << name << " " << saturation_index << " "
       << (*molality)[chemical_system_id] << " " << reaction_irreversibility
       << "\n";
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
