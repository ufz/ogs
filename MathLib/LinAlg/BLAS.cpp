/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BLAS.h"


#ifdef USE_PETSC

#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"

// Global PETScMatrix/PETScVector //////////////////////////////////////////

namespace MathLib { namespace BLAS
{

// Vector

void copy(PETScVector const& x, PETScVector& y)
{
    y = x;
}

void scale(PETScVector& x, double const a)
{
    VecScale(x.getRawVector(), a);
}

// y = a*y + X
void aypx(PETScVector& y, double const a, PETScVector const& x)
{
    // TODO check sizes
    VecAYPX(y.getRawVector(), a, x.getRawVector());
}

// y = a*x + y
void axpy(PETScVector& y, double const a, PETScVector const& x)
{
    // TODO check sizes
    VecAXPY(y.getRawVector(), a, x.getRawVector());
}

// y = a*x + y
void axpby(PETScVector& y, double const a, double const b, PETScVector const& x)
{
    // TODO check sizes
    VecAXPBY(y.getRawVector(), a, b, x.getRawVector());
}


// Matrix

void copy(PETScMatrix const& A, PETScMatrix& B)
{
    B = A;
}

// A = a*A
void scale(PETScMatrix& A, double const a)
{
    MatScale(A.getRawMatrix(), a);
}

// Y = a*Y + X
void aypx(PETScMatrix& Y, double const a, PETScMatrix const& X)
{
    // TODO check sizes
    // TODO sparsity pattern, currently they are assumed to be different (slow)
    MatAYPX(Y.getRawMatrix(), a, X.getRawMatrix(),
            DIFFERENT_NONZERO_PATTERN);
}

// Y = a*X + Y
void axpy(PETScMatrix& Y, double const a, PETScMatrix const& X)
{
    // TODO check sizes
    // TODO sparsity pattern, currently they are assumed to be different (slow)
    MatAXPY(Y.getRawMatrix(), a, X.getRawMatrix(),
            DIFFERENT_NONZERO_PATTERN);
}


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(PETScMatrix const& A, PETScVector const& x, PETScVector& y)
{
    // TODO check sizes
    assert(&x != &y);
    MatMult(A.getRawMatrix(), x.getRawVector(), y.getRawVector());
}

// v3 = A*v1 + v2
void matMultAdd(PETScMatrix const& A, PETScVector const& v1,
                       PETScVector const& v2, PETScVector& v3)
{
    // TODO check sizes
    assert(&v1 != &v3);
    MatMultAdd(A.getRawMatrix(), v1.getRawVector(), v2.getRawVector(), v3.getRawVector());
}

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

void copy(MEV const& x, MEV& y)
{
    y = x;
}

void scale(MEV& x, double const a)
{
    x *= a;
}

// y = a*y + X
void aypx(MEV& y, double const a, MEV const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a*y.getRawVector() + x.getRawVector();
}

// y = a*x + y
void axpy(MEV& y, double const a, MEV const& x)
{
    // TODO: does that break anything?
    y.getRawVector() += a*x.getRawVector();
}

// y = a*x + y
void axpby(MEV& y, double const a, double const b, MEV const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a*x.getRawVector() + b*y.getRawVector();
}


// Matrix

void copy(MEM const& A, MEM& B)
{
    B = A;
}

// A = a*A
void scale(MEM& A, double const a)
{
    // TODO: does that break anything?
    A.getRawMatrix() *= a;
}

// Y = a*Y + X
void aypx(MEM& Y, double const a, MEM const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a*Y.getRawMatrix() + X.getRawMatrix();
}

// Y = a*X + Y
void axpy(MEM& Y, double const a, MEM const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a*X.getRawMatrix() + Y.getRawMatrix();
}


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(MEM const& A, MEV const& x, MEV& y)
{
    assert(&x != &y);
    A.multiply(x, y);
}

// v3 = A*v1 + v2
void matMultAdd(MEM const& A, MEV const& v1, MEV const& v2, MEV& v3)
{
    assert(&v1 != &v3);
    // TODO: does that break anything?
    v3.getRawVector() = v2.getRawVector() + A.getRawMatrix()*v1.getRawVector();
}

} // namespace BLAS

} // namespace MathLib

#endif
