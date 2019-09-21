/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iosfwd>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Knobs
{
    Knobs(int max_iterations_, double relative_convergence_tolerance_,
          double tolerance_, int step_size_, bool scaling_)
        : max_iterations(max_iterations_),
          relative_convergence_tolerance(relative_convergence_tolerance_),
          tolerance(tolerance_),
          step_size(step_size_),
          scaling(scaling_)
    {
    }

    friend std::ostream& operator<<(std::ostream& os, Knobs const& knobs);

    int const max_iterations;
    double const relative_convergence_tolerance;
    double const tolerance;
    int const step_size;
    bool const scaling;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
