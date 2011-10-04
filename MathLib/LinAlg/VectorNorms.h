/*
 * VectorNorms.h
 *
 *  Created on: Jun 6, 2011
 *      Author: TF
 */

#ifndef VECTORNORMS_H_
#define VECTORNORMS_H_

#include <cmath>

#include "MathTools.h"

namespace MathLib {

double normEuklid (double const * const vec, size_t n)
{
	return sqrt (scpr (vec, vec, n));
}

} // end namespace MathLib

#endif /* VECTORNORMS_H_ */
