// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

    Set_initial_data(&initial_aqueous_solution);
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
