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

#ifdef OGS_USE_EIGEN

// Dense Eigen matrix/vector //////////////////////////////////////////

#include <Eigen/Core>

namespace MathLib
{

void setMatrix(Eigen::MatrixXd& m,
               Eigen::MatrixXd::Index const rows, Eigen::MatrixXd::Index const cols,
               std::initializer_list<double> values);

void setMatrix(Eigen::MatrixXd& m, Eigen::MatrixXd const& tmp);

void addToMatrix(Eigen::MatrixXd& m,
                 Eigen::MatrixXd::Index const rows, Eigen::MatrixXd::Index const cols,
                 std::initializer_list<double> values);

double norm(Eigen::VectorXd const& x);

void setVector(Eigen::VectorXd& v, std::initializer_list<double> values);

} // namespace MathLib

#endif // OGS_USE_EIGEN


#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////

#include "MathLib/LinAlg/PETSc/PETScMatrix.h"

namespace MathLib
{

class PETScVector;

double norm(PETScVector const& x);

void setVector(PETScVector& v,
                      std::initializer_list<double> values);


void setMatrix(PETScMatrix& m,
               PETScMatrix::IndexType const rows,
               PETScMatrix::IndexType const cols,
               std::initializer_list<double> values);

void setMatrix(PETScMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(PETScMatrix& m,
                 PETScMatrix::IndexType const rows,
                 PETScMatrix::IndexType const cols,
                 std::initializer_list<double> values);

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////

#include "MathLib/LinAlg/Eigen/EigenMatrix.h"

namespace MathLib
{

class EigenVector;

void setVector(EigenVector& v,
                      std::initializer_list<double> values);

void setMatrix(EigenMatrix& m,
               EigenMatrix::IndexType const rows,
               EigenMatrix::IndexType const cols,
               std::initializer_list<double> values);

void setMatrix(EigenMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(EigenMatrix& m,
                 EigenMatrix::IndexType const rows,
                 EigenMatrix::IndexType const cols,
                 std::initializer_list<double> values);

} // namespace MathLib

#endif // OGS_USE_EIGEN

#endif // MATHLIB_UNIFIED_MATRIX_SETTERS_H
