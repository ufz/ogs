/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ComputeResiduum.h"

#include "MathLib/LinAlg/LinAlg.h"

namespace ProcessLib
{
GlobalVector computeResiduum(GlobalVector const& x, GlobalVector const& xdot,
                             GlobalMatrix const& M, GlobalMatrix const& K,
                             GlobalVector const& b)
{
    using namespace MathLib::LinAlg;
    GlobalVector residuum;
    matMult(M, xdot, residuum);
    matMultAdd(K, x, residuum, residuum);
    axpy(residuum, -1, b);
    scale(residuum, -1);
    return residuum;
}
}  // namespace ProcessLib
