// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <iosfwd>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Knobs
{
    friend std::ostream& operator<<(std::ostream& os, Knobs const& knobs);

    int const max_iterations;
    double const relative_convergence_tolerance;
    double const tolerance;
    int const step_size;
    bool const scaling;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
