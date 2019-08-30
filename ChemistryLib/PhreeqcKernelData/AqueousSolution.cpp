/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AqueousSolution.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
AqueousSolution::AqueousSolution(double const temperature,
                                 double const pressure, double const pe_value,
                                 cxxISolution const& initial_aqueous_solution)
{
    new_def = true;
    tc = temperature;
    patm = pressure;
    pe = pe_value;
    initial_data = new cxxISolution(initial_aqueous_solution);
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
