/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "BaseLib/Error.h"
#include "LinAlgEnums.h"

#ifdef USE_PETSC
#include <petscsystypes.h>
#endif

namespace MathLib
{

/*! \namespace MathLib::LinAlg
 * Some general linear algebra functionality.
 *
 * By using the provided functions linear algebra capabilities can be
 * used for different matrix and vector types in a way that is agnostic
 * towards the specific type used.
 *
 * For documentation, refer to that of the templated method. All specializations
 * or overload must behave in the same way.
 */
namespace LinAlg
{

// Matrix or Vector

//! Copies \c x to \c y.
template<typename MatrixOrVector>
void copy(MatrixOrVector const& x, MatrixOrVector& y)
{
    y = x;
}

//! Scales \c x by \c a
template <typename MatrixOrVector>
void scale(MatrixOrVector& x, double const a)
{
    x *= a;
}

//! Computes \f$ y = a \cdot y + x \f$.
template <typename MatrixOrVector>
void aypx(MatrixOrVector& y, double const a, MatrixOrVector const& x)
{
    y = a*y + x;
}

//! Computes \f$ y = a \cdot x + y \f$.
template <typename MatrixOrVector>
void axpy(MatrixOrVector& y, double const a, MatrixOrVector const& x)
{
    y += a*x;
}

//! Computes \f$ y = a \cdot x + b \cdot y \f$.
template <typename MatrixOrVector>
void axpby(MatrixOrVector& y, double const a, double const b,
           MatrixOrVector const& x)
{
    y = a*x + b*y;
}

//! Computes \f$w = x/y\f$ componentwise.
template<typename MatrixOrVector>
void componentwiseDivide(MatrixOrVector& w,
                         MatrixOrVector const& x, MatrixOrVector const& y);

//! Computes the Manhattan norm of \c x.
template<typename MatrixOrVector>
double norm1(MatrixOrVector const& x);

//! Computes the Euclidean norm of \c x.
template<typename MatrixOrVector>
double norm2(MatrixOrVector const& x);

//! Computes the maximum norm of \c x.
template<typename MatrixOrVector>
double normMax(MatrixOrVector const& x);

template<typename MatrixOrVector>
double norm(MatrixOrVector const& x, MathLib::VecNormType type)
{
    switch (type) {
        case MathLib::VecNormType::NORM1:
            return norm1(x);
        case MathLib::VecNormType::NORM2:
            return norm2(x);
        case MathLib::VecNormType::INFINITY_N:
            return normMax(x);
        default:
            OGS_FATAL("Invalid norm type given.");
    }
}

template <typename Matrix>
void finalizeAssembly(Matrix& /*A*/);

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

}  // namespace LinAlg
}  // namespace MathLib

// Global PETScMatrix/PETScVector //////////////////////////////////////////
#ifdef USE_PETSC

namespace MathLib {

class PETScMatrix;
class PETScVector;

namespace LinAlg
{

// Vector

/// Set local accessible vector in order to get entries.
/// Call this before call operator[] or get(...) of x.
void setLocalAccessibleVector(PETScVector const& x);

void set(PETScVector& x, PetscScalar const a);

void copy(PETScVector const& x, PETScVector& y);

void scale(PETScVector& x, PetscScalar const a);

// y = a*y + X
void aypx(PETScVector& y, PetscScalar const a, PETScVector const& x);

// y = a*x + y
void axpy(PETScVector& y, PetscScalar const a, PETScVector const& x);

// y = a*x + b*y
void axpby(PETScVector& y, PetscScalar const a, PetscScalar const b,
           PETScVector const& x);

// Matrix

void copy(PETScMatrix const& A, PETScMatrix& B);

// A = a*A
void scale(PETScMatrix& A, PetscScalar const a);

// Y = a*Y + X
void aypx(PETScMatrix& Y, PetscScalar const a, PETScMatrix const& X);

// Y = a*X + Y
void axpy(PETScMatrix& Y, PetscScalar const a, PETScMatrix const& X);

// Matrix and Vector

// v3 = A*v1 + v2
void matMult(PETScMatrix const& A, PETScVector const& x, PETScVector& y);

// y = A*x
void matMultAdd(PETScMatrix const& A, PETScVector const& v1,
                       PETScVector const& v2, PETScVector& v3);

void finalizeAssembly(PETScMatrix& A);
void finalizeAssembly(PETScVector& x);

}} // namespaces


// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////
#else

namespace MathLib {

class EigenMatrix;
class EigenVector;

namespace LinAlg
{

// Vector

/**
    Set local accessible vector in order to get entries. Call this
    before call operator[] or get(...) of x.
    The function only has computation if DDC is enabled,
    e.g. parallel computing. Up to now, Eigen vector is not used for
    global vectors in parallel computing. Therefore this function is
    empty in the current status.
*/
void setLocalAccessibleVector(EigenVector const& x);

void set(EigenVector& x, double const a);

void copy(EigenVector const& x, EigenVector& y);

void scale(EigenVector& x, double const a);

// y = a*y + X
void aypx(EigenVector& y, double const a, EigenVector const& x);

// y = a*x + y
void axpy(EigenVector& y, double const a, EigenVector const& x);

// y = a*x + b*y
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

void finalizeAssembly(EigenMatrix& x);
void finalizeAssembly(EigenVector& A);

} // namespace LinAlg

} // namespace MathLib

#endif

namespace MathLib::LinAlg
{

/**
 * Compute the relative norm between two vectors by \f$ \|x-y\|/\|x\| \f$.
 *
 * \param x Vector x
 * \param y Vector y
 * \param norm_type The norm type of global vector
 * \return \f$ \|x-y\|/\|x\| \f$.
 */
template <typename VectorType>
double computeRelativeNorm(VectorType const& x, VectorType const& y,
                           MathLib::VecNormType norm_type)
{
    if (norm_type == MathLib::VecNormType::INVALID)
    {
        OGS_FATAL("An invalid norm type has been passed");
    }

    // Stores \f$diff = x-y\f$.
    VectorType diff;
    MathLib::LinAlg::copy(x, diff);
    MathLib::LinAlg::axpy(diff, -1.0, y);

    const double norm_diff = MathLib::LinAlg::norm(diff, norm_type);

    const double norm_x = MathLib::LinAlg::norm(x, norm_type);
    if (norm_x > std::numeric_limits<double>::epsilon())
    {
        return norm_diff / norm_x;
    }

    // Both of norm_x and norm_diff are close to zero
    if (norm_diff < std::numeric_limits<double>::epsilon())
    {
        return 1.0;
    }

    // Only norm_x is close to zero
    return norm_diff / std::numeric_limits<double>::epsilon();
}
}  // namespace MathLib::LinAlg
