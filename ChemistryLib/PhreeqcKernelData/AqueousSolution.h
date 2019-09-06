/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/common/phrqtype.h"
#include "ThirdParty/iphreeqc/src/src/phreeqcpp/Solution.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class AqueousSolution final : public cxxSolution
{
public:
    AqueousSolution(double const temperature, double const pressure,
                    double const pe_value,
                    cxxISolution const& initial_aqueous_solution);
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
