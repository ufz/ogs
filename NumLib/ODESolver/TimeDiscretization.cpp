/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on March 31, 2017, 10:30 AM
 */

#include "TimeDiscretization.h"

#include "MathLib/LinAlg/MatrixVectorTraits.h"

namespace NumLib
{
double computeRelativeChangeFromPreviousTimestep(GlobalVector const& x,
                                                 GlobalVector const& x_old,
                                                 MathLib::VecNormType norm_type)
{
    if (norm_type == MathLib::VecNormType::INVALID)
    {
        OGS_FATAL("An invalid norm type has been passed");
    }

    // Stores \f$ u_{n+1}-u_{n}\f$.
    GlobalVector dx;
    MathLib::LinAlg::copy(x, dx);  // copy x to dx.

    // dx = x - x_old --> x - dx --> dx
    MathLib::LinAlg::axpy(dx, -1.0, x_old);
    const double norm_dx = MathLib::LinAlg::norm(dx, norm_type);

    const double norm_x = MathLib::LinAlg::norm(x, norm_type);
    if (norm_x > std::numeric_limits<double>::epsilon())
    {
        return norm_dx / norm_x;
    }

    // Both of norm_x and norm_dx are close to zero
    if (norm_dx < std::numeric_limits<double>::epsilon())
    {
        return 1.0;
    }

    // Only norm_x is close to zero
    return norm_dx / std::numeric_limits<double>::epsilon();
}
}  // end of namespace NumLib
