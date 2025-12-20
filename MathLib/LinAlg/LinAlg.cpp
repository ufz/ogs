// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "LinAlg.h"

// TODO reorder LinAlg function signatures?

// Global PETScMatrix/PETScVector //////////////////////////////////////////
#ifdef USE_PETSC

#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"

namespace MathLib
{
namespace LinAlg
{
// Vector

void setLocalAccessibleVector(PETScVector const& x)
{
    x.setLocalAccessibleVector();
}

void set(PETScVector& x, PetscScalar const a)
{
    VecSet(x.getRawVector(), a);
}

void copy(PETScVector const& x, PETScVector& y)
{
    if (!y.getRawVector())
        y.shallowCopy(x);
    VecCopy(x.getRawVector(), y.getRawVector());
}

void scale(PETScVector& x, PetscScalar const a)
{
    VecScale(x.getRawVector(), a);
}

// y = a*y + X
void aypx(PETScVector& y, PetscScalar const a, PETScVector const& x)
{
    // TODO check sizes
    VecAYPX(y.getRawVector(), a, x.getRawVector());
}

// y = a*x + y
void axpy(PETScVector& y, PetscScalar const a, PETScVector const& x)
{
    // TODO check sizes
    VecAXPY(y.getRawVector(), a, x.getRawVector());
}

// y = a*x + b*y
void axpby(PETScVector& y, PetscScalar const a, PetscScalar const b,
           PETScVector const& x)
{
    // TODO check sizes
    VecAXPBY(y.getRawVector(), a, b, x.getRawVector());
}

// Explicit specialization
// Computes w = x/y componentwise.
// \note  that VecPointwiseDivide avoids to divide by values that are
// identically zero such as
//   for (int i=0; i<n; i++)
//   {
//      w[i] = y[i] == 0.0 ? 0.0 : x[i] / y[i];
//   }
//
template <>
void componentwiseDivide(PETScVector& w, PETScVector const& x,
                         PETScVector const& y)
{
    VecPointwiseDivide(w.getRawVector(), x.getRawVector(), y.getRawVector());
}

// Explicit specialization
// Computes the Manhattan norm of x
template <>
double norm1(PETScVector const& x)
{
    PetscScalar norm = 0.;
    VecNorm(x.getRawVector(), NORM_1, &norm);
    return norm;
}

// Explicit specialization
// Computes the Euclidean norm of x
template <>
double norm2(PETScVector const& x)
{
    PetscScalar norm = 0.;
    VecNorm(x.getRawVector(), NORM_2, &norm);
    return norm;
}

// Explicit specialization
// Computes the Maximum norm of x
template <>
double normMax(PETScVector const& x)
{
    PetscScalar norm = 0.;
    VecNorm(x.getRawVector(), NORM_INFINITY, &norm);
    return norm;
}

// Matrix

void copy(PETScMatrix const& A, PETScMatrix& B)
{
    B = A;
}

// A = a*A
void scale(PETScMatrix& A, PetscScalar const a)
{
    MatScale(A.getRawMatrix(), a);
}

// Y = a*Y + X
void aypx(PETScMatrix& Y, PetscScalar const a, PETScMatrix const& X)
{
    // TODO check sizes
    // TODO sparsity pattern, currently they are assumed to be different (slow)
    MatAYPX(Y.getRawMatrix(), a, X.getRawMatrix(), DIFFERENT_NONZERO_PATTERN);
}

// Y = a*X + Y
void axpy(PETScMatrix& Y, PetscScalar const a, PETScMatrix const& X)
{
    // TODO check sizes
    // TODO sparsity pattern, currently they are assumed to be different (slow)
    MatAXPY(Y.getRawMatrix(), a, X.getRawMatrix(), DIFFERENT_NONZERO_PATTERN);
}

// Matrix and Vector

// v3 = A*v1 + v2
void matMult(PETScMatrix const& A, PETScVector const& x, PETScVector& y)
{
    // TODO check sizes
    assert(&x != &y);
    if (!y.getRawVector())
        y.shallowCopy(x);
    MatMult(A.getRawMatrix(), x.getRawVector(), y.getRawVector());
}

// v3 = A*v1 + v2
void matMultAdd(PETScMatrix const& A, PETScVector const& v1,
                PETScVector const& v2, PETScVector& v3)
{
    // TODO check sizes
    assert(&v1 != &v3);
    if (!v3.getRawVector())
        v3.shallowCopy(v1);
    MatMultAdd(A.getRawMatrix(), v1.getRawVector(), v2.getRawVector(),
               v3.getRawVector());
}

void linearSysNormalize(PETScMatrix const& /*A*/, PETScMatrix& /*new_A*/,
                        PETScVector const& /*b*/, PETScVector& /*new_b*/)
{
    // The following block is deactivated, because there is no tests yet for the
    // normalization operation in PETSc. This will be a task for later.
    /*
    assert(&A != &new_A);
    assert(&b != &new_b);

    PetscInt n_rows(0);
    PetscInt n_cols(0);
    MatGetSize(A.getRawMatrix(), &n_rows, &n_cols);
    // only when A matrix is not square
    if (n_rows != n_cols)
    {
        // new_b = A^T * b
        MatMultTranspose(A.getRawMatrix(), b.getRawVector(),
                         new_b.getRawVector());
        // new_A = A^T * A
        MatTranspose(A.getRawMatrix(), MAT_INITIAL_MATRIX,
                     &(new_A.getRawMatrix()));
    }
    */
    OGS_FATAL(
        "Normalization operation is not implemented yet for PETSc library! "
        "Program terminated.");
}

void finalizeAssembly(PETScMatrix& A)
{
    A.finalizeAssembly(MAT_FINAL_ASSEMBLY);
}

void finalizeAssembly(PETScVector& x)
{
    x.finalizeAssembly();
}

}  // namespace LinAlg
}  // namespace MathLib

// Sparse global EigenMatrix/EigenVector
// //////////////////////////////////////////
#else

#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"

namespace MathLib
{
namespace LinAlg
{
// Vector

void setLocalAccessibleVector(EigenVector const& /*x*/) {}

void set(EigenVector& x, double const a)
{
    x.getRawVector().setConstant(a);
}

void copy(EigenVector const& x, EigenVector& y)
{
    y = x;
}

void scale(EigenVector& x, double const a)
{
    x.getRawVector() *= a;
}

// y = a*y + X
void aypx(EigenVector& y, double const a, EigenVector const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a * y.getRawVector() + x.getRawVector();
}

// y = a*x + y
void axpy(EigenVector& y, double const a, EigenVector const& x)
{
    // TODO: does that break anything?
    y.getRawVector() += a * x.getRawVector();
}

// y = a*x + b*y
void axpby(EigenVector& y, double const a, double const b, EigenVector const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a * x.getRawVector() + b * y.getRawVector();
}

// Explicit specialization
// Computes w = x/y componentwise.
template <>
void componentwiseDivide(EigenVector& w, EigenVector const& x,
                         EigenVector const& y)
{
    w.getRawVector().noalias() = x.getRawVector().binaryExpr(
        y.getRawVector(),
        [](auto const x, auto const y) { return y == 0 ? 0.0 : x / y; });
}

// Explicit specialization
// Computes the Manhattan norm of x
template <>
double norm1(EigenVector const& x)
{
    return x.getRawVector().lpNorm<1>();
}

// Explicit specialization
// Euclidean norm
template <>
double norm2(EigenVector const& x)
{
    return x.getRawVector().norm();
}

// Explicit specialization
// Computes the Maximum norm of x
template <>
double normMax(EigenVector const& x)
{
    return x.getRawVector().lpNorm<Eigen::Infinity>();
}

// Matrix

void copy(EigenMatrix const& A, EigenMatrix& B)
{
    B = A;
}

// A = a*A
void scale(EigenMatrix& A, double const a)
{
    // TODO: does that break anything?
    A.getRawMatrix() *= a;
}

// Y = a*Y + X
void aypx(EigenMatrix& Y, double const a, EigenMatrix const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a * Y.getRawMatrix() + X.getRawMatrix();
}

// Y = a*X + Y
void axpy(EigenMatrix& Y, double const a, EigenMatrix const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a * X.getRawMatrix() + Y.getRawMatrix();
}

// Matrix and Vector

// v3 = A*v1 + v2
void matMult(EigenMatrix const& A, EigenVector const& x, EigenVector& y)
{
    assert(&x != &y);
    y.getRawVector() = A.getRawMatrix() * x.getRawVector();
}

// v3 = A*v1 + v2
void matMultAdd(EigenMatrix const& A, EigenVector const& v1,
                EigenVector const& v2, EigenVector& v3)
{
    assert(&v1 != &v3);
    // TODO: does that break anything?
    v3.getRawVector() =
        v2.getRawVector() + A.getRawMatrix() * v1.getRawVector();
}

void linearSysNormalize(EigenMatrix const& A, EigenMatrix& new_A,
                        EigenVector const& b, EigenVector& new_b)
{
    // make sure that new_A and new_b are not the same memory
    assert(&A != &new_A);
    assert(&b != &new_b);

    if (A.getRawMatrix().rows() == A.getRawMatrix().cols())
    {
        WARN(
            "The number of rows and columns are the same for the LHS matrix."
            "Are you sure you still need to normalize the LHS matrix and RHS "
            "vector? ");
    }

    new_b.getRawVector() = A.getRawMatrix().transpose() * b.getRawVector();
    new_A.getRawMatrix() = A.getRawMatrix().transpose() * A.getRawMatrix();
}

void finalizeAssembly(EigenMatrix& x)
{
    x.getRawMatrix().makeCompressed();
}

void finalizeAssembly(EigenVector& /*x*/) {}

}  // namespace LinAlg

}  // namespace MathLib

#endif
