/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/Error.h"
#include "EigenMatrix.h" // for EigenMatrix::IndexType

namespace MathLib
{
class EigenVector;

/**
 * apply known solutions to a system of linear equations
 *
 * @param A                 Coefficient matrix
 * @param b                 RHS vector
 * @param vec_knownX_id_    a vector of known solution entry IDs
 * @param vec_knownX_x_     a vector of known solutions
 * @param penalty_scaling value for scaling some matrix and right hand side
 * entries to enforce some conditions
 */
void applyKnownSolution(EigenMatrix &A, EigenVector &b, EigenVector &/*x*/,
        const std::vector<EigenMatrix::IndexType> &vec_knownX_id_,
        const std::vector<double> &vec_knownX_x_, double penalty_scaling = 1e+10);

inline
void applyKnownSolution(Eigen::MatrixXd const &A, Eigen::VectorXd &b, Eigen::VectorXd &/*x*/,
        const std::vector<Eigen::MatrixXd::Index> &vec_knownX_id_,
        const std::vector<double> &vec_knownX_x_, double penalty_scaling = 1e+10)
{
    (void) A; (void) b; (void) vec_knownX_id_; (void) vec_knownX_x_;
    (void) penalty_scaling;

    OGS_FATAL("Method not implemented."); // TODO implement
}

}  // namespace MathLib
