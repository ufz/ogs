/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 *
 * \date   2011-06-06 -- 2013-12-10
 *  
 * \brief  Definition of vector norm functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTORNORMS_H_
#define VECTORNORMS_H_

#include <cmath>

#include "MathTools.h"

namespace MathLib {

/// Norm type. Not declared as class type in order to use the members as integers. 
enum VectorNormType 
{
   sum_abs_entries, ///< \f$\sum_i |x_i|\f$
   euclidean,       ///< \f$\sqrt(\sum_i (x_i)^2)\f$
   max_abs_entry    ///< \f$\mathrm{max}_i |x_i|\f$
};

inline double normEuklid (double const * const vec, std::size_t n)
{
	return sqrt (scpr (vec, vec, n));
}

} // end namespace MathLib

#endif /* VECTORNORMS_H_ */
