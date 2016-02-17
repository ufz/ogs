/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_BLAS_H
#define MATHLIB_BLAS_H

#ifdef OGS_USE_EIGEN

#include <Eigen/Core>

namespace MathLib { namespace BLAS
{

// Dense Eigen matrix/vector //////////////////////////////////////////
// TODO change to templates

using EM = Eigen::MatrixXd;
using EV = Eigen::VectorXd;


// Vector

void copy(EV const& x, EV& y);

void scale(EV& x, double const a);

// y = a*y + X
void aypx(EV& y, double const a, EV const& x);

// y = a*x + y
void axpy(EV& y, double const a, EV const& x);

// y = a*x + y
void axpby(EV& y, double const a, double const b, EV const& x);


// Matrix

void copy(EM const& A, EM& B);

// A = a*A
void scale(EM& A, double const a);

// Y = a*Y + X
void aypx(EM& Y, double const a, EM const& X);

// Y = a*X + Y
void axpy(EM& Y, double const a, EM const& X);


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(EM const& A, EV const& x, EV& y);

// v3 = A*v1 + v2
void matMultAdd(EM const& A, EV const& v1, EV const& v2, EV& v3);

}} // namespaces

#endif


#ifdef USE_PETSC

#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"

// Global PETScMatrix/PETScVector //////////////////////////////////////////

namespace MathLib { namespace BLAS
{

// Vector

void copy(PETScVector const& x, PETScVector& y);

void scale(PETScVector& x, double const a);

// y = a*y + X
void aypx(PETScVector& y, double const a, PETScVector const& x);

// y = a*x + y
void axpy(PETScVector& y, double const a, PETScVector const& x);

// y = a*x + y
void axpby(PETScVector& y, double const a, double const b, PETScVector const& x);


// Matrix

void copy(PETScMatrix const& A, PETScMatrix& B);

// A = a*A
void scale(PETScMatrix& A, double const a);

// Y = a*Y + X
void aypx(PETScMatrix& Y, double const a, PETScMatrix const& X);

// Y = a*X + Y
void axpy(PETScMatrix& Y, double const a, PETScMatrix const& X);


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(PETScMatrix const& A, PETScVector const& x, PETScVector& y);

// v3 = A*v1 + v2
void matMultAdd(PETScMatrix const& A, PETScVector const& v1,
                       PETScVector const& v2, PETScVector& v3);

}} // namespaces


#elif defined(OGS_USE_EIGEN)

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////

#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"


namespace MathLib { namespace BLAS
{

using MEM = MathLib::EigenMatrix;
using MEV = MathLib::EigenVector;


// Vector

void copy(MEV const& x, MEV& y);

void scale(MEV& x, double const a);

// y = a*y + X
void aypx(MEV& y, double const a, MEV const& x);

// y = a*x + y
void axpy(MEV& y, double const a, MEV const& x);

// y = a*x + y
void axpby(MEV& y, double const a, double const b, MEV const& x);


// Matrix

void copy(MEM const& A, MEM& B);

// A = a*A
void scale(MEM& A, double const a);

// Y = a*Y + X
void aypx(MEM& Y, double const a, MEM const& X);

// Y = a*X + Y
void axpy(MEM& Y, double const a, MEM const& X);


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(MEM const& A, MEV const& x, MEV& y);

// v3 = A*v1 + v2
void matMultAdd(MEM const& A, MEV const& v1, MEV const& v2, MEV& v3);

} // namespace BLAS

} // namespace MathLib

#endif

#endif // MATHLIB_BLAS_H
