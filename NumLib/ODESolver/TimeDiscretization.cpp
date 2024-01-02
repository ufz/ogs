/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 31, 2017, 10:30 AM
 */

#include "TimeDiscretization.h"

#include "MathLib/LinAlg/MatrixVectorTraits.h"

namespace NumLib
{
void BackwardEuler::getWeightedOldX(GlobalVector& y,
                                    GlobalVector const& x_old) const
{
    namespace LinAlg = MathLib::LinAlg;

    // y = x_old / delta_t
    LinAlg::copy(x_old, y);
    LinAlg::scale(y, 1.0 / _delta_t);
}

}  // end of namespace NumLib
