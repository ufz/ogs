/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file VectorNorms.h
 *
 * Created on 2011-06-06 by Thomas Fischer
 */

#ifndef VECTORNORMS_H_
#define VECTORNORMS_H_

#include <cmath>

#include "MathTools.h"

namespace MathLib {

double normEuklid (double const * const vec, std::size_t n)
{
	return sqrt (scpr (vec, vec, n));
}

} // end namespace MathLib

#endif /* VECTORNORMS_H_ */
