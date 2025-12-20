// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Knobs.h"

#include <ostream>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::ostream& operator<<(std::ostream& os, Knobs const& knobs)
{
    os << "KNOBS"
       << "\n";
    os << "-iterations " << knobs.max_iterations << "\n";
    os << "-convergence_tolerance " << knobs.relative_convergence_tolerance
       << "\n";
    os << "-tolerance " << knobs.tolerance << "\n";
    os << "-step_size " << knobs.step_size << "\n";
    os << "-diagonal_scale " << knobs.scaling << "\n";

    return os;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
