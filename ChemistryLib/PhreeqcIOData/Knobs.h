/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    friend std::ostream& operator<<(std::ostream& os, Knobs const& knobs);

    int const max_iterations;
    double const relative_convergence_tolerance;
    double const tolerance;
    int const step_size;
    bool const scaling;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
