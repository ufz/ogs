/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TESOGS5MaterialModels.h"

namespace ProcessLib
{
namespace TES
{
const double FluidHeatConductivityN2::A[5] = {0.46649, -0.57015, 0.19164,
                                              -0.03708, 0.00241};

const double FluidHeatConductivityN2::f[9] = {
    -0.837079888737e3,   0.37914711487e2,    -0.601737844275,
    0.350418363823e1,    -0.874955653028e-5, 0.148968607239e-7,
    -0.256370354277e-11, 0.100773735767e1,   0.335340610e4};

const double FluidHeatConductivityN2::C[4] = {3.3373542, 0.37098251, 0.89913456,
                                              0.16972505};

const double FluidViscosityN2::A[5] = {0.46649, -0.57015, 0.19164, -0.03708,
                                       0.00241};
const double FluidViscosityN2::C[5] = {-20.09997, 3.4376416, -1.4470051,
                                       -0.027766561, -0.21662362};

const double FluidHeatConductivityH2O::a[4] = {0.0102811, 0.0299621, 0.0156146,
                                               -0.00422464};

}  // TES
}  // ProcessLib
