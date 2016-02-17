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

#include<cassert>


namespace MathLib
{

/*! \namespace MathLib::BLAS
 * Some general linear algebra functionality.
 *
 * By using the provided functions linear algebra capabilities can be
 * used for different matrix and vector types in a way that is agnostic
 * towards the specific type used.
 *
 * For documentation, refer to that of the templated method. All specializations
 * or overload must behave in the same way.
 */
namespace BLAS
{

// Matrix or Vector

//! Copies \c x to \c y.
template<typename MatrixOrVector>
void copy(MatrixOrVector const& x, MatrixOrVector& y)
{
    y = x;
}

//! Scales \c x by \c a
template<typename MatrixOrVector>
void scale(MatrixOrVector& x, double const a)
{
    x *= a;
}

//! Computes \f$ y = a \cdot y + x \f$.
template<typename MatrixOrVector>
void aypx(MatrixOrVector& y, double const a, MatrixOrVector const& x)
{
    y = a*y + x;
}

//! Computes \f$ y = a \cdot x + y \f$.
template<typename MatrixOrVector>
void axpy(MatrixOrVector& y, double const a, MatrixOrVector const& x)
{
    y += a*x;
}

//! Computes \f$ y = a \cdot x + b \cdot y \f$.
template<typename MatrixOrVector>
void axpby(MatrixOrVector& y, double const a, double const b, MatrixOrVector const& x)
{
    y = a*x + b*y;
}


// Matrix and Vector

/*! Computes \f$ y = A \cdot x \f$.
 *
 * \note \c x must not be the same object as \c y.
 *       This restirction has been chosen in order to fulfill
 *       the requirements of the respective PETSc function.
 */
template<typename Matrix, typename Vector>
void matMult(Matrix const& A, Vector const& x, Vector& y)
{
    assert(&x != &y);
    y = A*x;
}

/*! Computes \f$ v_3 = A \cdot v_1 + v_2 \f$.
 *
 * \note \c x must not be the same object as \c y.
 *       This restirction has been chosen in order to fulfill
 *       the requirements of the respective PETSc function.
 */
template<typename Matrix, typename Vector>
void matMultAdd(Matrix const& A, Vector const& v1, Vector const& v2, Vector& v3)
{
    assert(&v1 != &v3);
    v3 = v2 + A*v1;
}

}} // namespaces


// Global PETScMatrix/PETScVector //////////////////////////////////////////
#ifdef USE_PETSC

namespace MathLib {

class PETScMatrix;
class PETScVector;

namespace BLAS
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

// y = A*x
void matMultAdd(PETScMatrix const& A, PETScVector const& v1,
                       PETScVector const& v2, PETScVector& v3);

}} // namespaces



// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////
#elif defined(OGS_USE_EIGEN)

namespace MathLib {

class EigenMatrix;
class EigenVector;

namespace BLAS
{

// Vector

void copy(EigenVector const& x, EigenVector& y);

void scale(EigenVector& x, double const a);

// y = a*y + X
void aypx(EigenVector& y, double const a, EigenVector const& x);

// y = a*x + y
void axpy(EigenVector& y, double const a, EigenVector const& x);

// y = a*x + y
void axpby(EigenVector& y, double const a, double const b, EigenVector const& x);


// Matrix

void copy(EigenMatrix const& A, EigenMatrix& B);

// A = a*A
void scale(EigenMatrix& A, double const a);

// Y = a*Y + X
void aypx(EigenMatrix& Y, double const a, EigenMatrix const& X);

// Y = a*X + Y
void axpy(EigenMatrix& Y, double const a, EigenMatrix const& X);


// Matrix and Vector

// y = A*x
void matMult(EigenMatrix const& A, EigenVector const& x, EigenVector& y);

// v3 = A*v1 + v2
void matMultAdd(EigenMatrix const& A, EigenVector const& v1,
                EigenVector const& v2, EigenVector& v3);

} // namespace BLAS

} // namespace MathLib

#endif

#endif // MATHLIB_BLAS_H
