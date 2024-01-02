/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ComputeResiduum.h"

#include "MathLib/LinAlg/LinAlg.h"

namespace ProcessLib
{
GlobalVector computeResiduum(double const dt, GlobalVector const& x,
                             GlobalVector const& x_prev, GlobalMatrix const& M,
                             GlobalMatrix const& K, GlobalVector const& b)
{
    using namespace MathLib::LinAlg;
    GlobalVector residuum;
    GlobalVector x_dot;
    copy(x, x_dot);                        // tmp = x
    axpy(x_dot, -1., x_prev);              // tmp = x - x_prev
    scale(x_dot, 1. / dt);                 // tmp = (x - x_prev)/dt
    matMult(M, x_dot, residuum);           // r = M*x_dot
    matMultAdd(K, x, residuum, residuum);  // r = M*x_dot + K*x
    axpy(residuum, -1., b);                // r = M*x_dot + K*x - b
    scale(residuum, -1.);                  // r = -r
    return residuum;
}
}  // namespace ProcessLib
