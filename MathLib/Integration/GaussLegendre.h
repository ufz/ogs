/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include "WeightedSum.h"

namespace MathLib
{

/// Gauss-Legendre quadrature method
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendre {
    static const unsigned Order = ORDER;
    static const double X[Order];
    static const double W[Order];
};

}   // namespace MathLib

#endif // GAUSSLEGENDRE_H_

