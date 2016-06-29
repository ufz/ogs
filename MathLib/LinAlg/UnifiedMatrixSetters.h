/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// TODO merge that with MatrixVectorTraits?

#ifndef MATHLIB_UNIFIED_MATRIX_SETTERS_H
#define MATHLIB_UNIFIED_MATRIX_SETTERS_H

#include <initializer_list>
#include "MatrixVectorTraits.h"

#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////
#include "PETSc/PETScVector.h"

namespace MathLib
{

class PETScMatrix;

double norm(PETScVector const& x);

void setVector(PETScVector& v,
               std::initializer_list<double> values);

void setVector(PETScVector& v, MatrixVectorTraits<PETScVector>::Index const index,
               double const value);

void setMatrix(PETScMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(PETScMatrix& m,
                 std::initializer_list<double> values);

void setMatrix(PETScMatrix& m,
               std::initializer_list<double> values);

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////
#include "Eigen/EigenVector.h"

namespace MathLib
{

class EigenMatrix;

void setVector(EigenVector& v,
               std::initializer_list<double> values);

void setVector(EigenVector& v, MatrixVectorTraits<EigenVector>::Index const index,
               double const value);

void setMatrix(EigenMatrix& m,
               std::initializer_list<double> values);

void setMatrix(EigenMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(EigenMatrix& m,
                 std::initializer_list<double> values);

} // namespace MathLib

#endif // OGS_USE_EIGEN

#endif // MATHLIB_UNIFIED_MATRIX_SETTERS_H
