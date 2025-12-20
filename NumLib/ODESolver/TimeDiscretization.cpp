// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
