/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of Lis utility functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISTOOLS_H_
#define LISTOOLS_H_

#include <vector>

namespace MathLib
{
class LisMatrix;
class LisVector;

/**
 * check Lis error codes
 *
 * @param err   Lis error code
 * @return success or not
 */
bool checkLisError(int err);

/**
 * apply known solutions to a system of linear equations
 *
 * This function introduces the constants into the system by the penalty method.
 *
 * @param A                 Coefficient matrix
 * @param b                 RHS vector
 * @param _vec_knownX_id    a vector of known solution entry IDs
 * @param _vec_knownX_x     a vector of known solutions
 */
void applyKnownSolution(LisMatrix &A, LisVector &b, const std::vector<std::size_t> &_vec_knownX_id, const std::vector<double> &_vec_knownX_x);

} // MathLib

#endif //LISTOOLS_H_

