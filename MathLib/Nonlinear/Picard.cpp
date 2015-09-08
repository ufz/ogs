/*
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Picard.h"

namespace MathLib
{
namespace Nonlinear
{
Picard::Picard()
    : _normType(VecNormType::INFINITY_N),
      _abs_tol(std::numeric_limits<double>::max()),
      _rel_tol(1e-6),
      _max_itr(25),
      _printErrors(false),
      _n_iterations(0),
      _abs_error(.0),
      _rel_error(.0)
{
}

}  // namespace Nonlinear

}  // namespace MathLib
