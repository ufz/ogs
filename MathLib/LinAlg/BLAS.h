#pragma once

#include<cassert>


#include <Eigen/Core>

namespace MathLib { namespace BLAS
{

// Dense Eigen matrix/vector //////////////////////////////////////////
// TODO change t otemplates

using EM = Eigen::MatrixXd;
using EV = Eigen::VectorXd;


// Vector

inline void copy(EV const& x, EV& y)
{
    y = x;
}

inline void scale(EV& x, double const a)
{
    x *= a;
}

// y = a*y + X
inline void aypx(EV& y, double const a, EV const& x)
{
    y = a*y + x;
}

// y = a*x + y
inline void axpy(EV& y, double const a, EV const& x)
{
    y += a*x;
}

// y = a*x + y
inline void axpby(EV& y, double const a, double const b, EV const& x)
{
    y = a*x + b*y;
}


// Matrix

inline void copy(EM const& A, EM& B)
{
    B = A;
}

// A = a*A
inline void scale(EM& A, double const a)
{
    A *= a;
}

// Y = a*Y + X
inline void aypx(EM& Y, double const a, EM const& X)
{
    Y = a*Y + X;
}

// Y = a*X + Y
inline void axpy(EM& Y, double const a, EM const& X)
{
    Y = a*X + Y;
}


// Matrix and Vector

// v3 = A*v1 + v2
inline void matMult(EM const& A, EV const& x, EV& y)
{
    assert(&x != &y);
    y = A*x;
}

// v3 = A*v1 + v2
inline void matMultAdd(EM const& A, EV const& v1, EV const& v2, EV& v3)
{
    assert(&v1 != &v3);
    v3 = v2 + A*v1;
}

}} // namespaces


#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////
///
///
namespace MathLib { namespace BLAS
{

// Vector

inline void copy(PETScVector const& x, PETScVector& y)
{
    y = x;
}

inline void scale(PETScVector& x, double const a)
{
    (void) x; (void) a;
    // TODO implement
}

// y = a*y + X
inline void aypx(PETScVector& y, double const a, PETScVector const& x)
{
    (void) y; (void) a; (void) x;
    // TODO implement
}

// y = a*x + y
inline void axpy(PETScVector& y, double const a, PETScVector const& x)
{
    (void) y; (void) a; (void) x;
    // TODO implement
}

// y = a*x + y
inline void axpby(PETScVector& y, double const a, double const b, PETScVector const& x)
{
    (void) y; (void) a; (void) b; (void) x;
    // TODO implement
}


// Matrix

inline void copy(PETScMatrix const& A, PETScMatrix& B)
{
    B = A;
}

// A = a*A
inline void scale(PETScMatrix& A, double const a)
{
    (void) A; (void) a;
    // TODO implement
}

// Y = a*Y + X
inline void aypx(PETScMatrix& Y, double const a, PETScMatrix const& X)
{
    (void) Y;(void) a; (void) X;
    // TODO implement
}

// Y = a*X + Y
inline void axpy(PETScMatrix& Y, double const a, PETScMatrix const& X)
{
    (void) Y; (void) a; (void) X;
    // TODO implement
}


// Matrix and Vector

// v3 = A*v1 + v2
inline void matMult(PETScMatrix const& A, PETScVector const& x, PETScVector& y)
{
    (void) A; (void) x; (void) y;
    assert(&x != &y);
    // TODO implement
}

// v3 = A*v1 + v2
inline void matMultAdd(PETScMatrix const& A, PETScVector const& v1,
                       PETScVector const& v2, PETScVector& v3)
{
    (void) A; (void) v1; (void) v2; (void) v3;
    assert(&v1 != &v3);
    // TODO implement
}

}} // namespaces


#else

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////

#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"


namespace MathLib { namespace BLAS
{

using MEM = MathLib::EigenMatrix;
using MEV = MathLib::EigenVector;


// Vector

inline void copy(MEV const& x, MEV& y)
{
    y = x;
}

inline void scale(MEV& x, double const a)
{
    x *= a;
}

// y = a*y + X
inline void aypx(MEV& y, double const a, MEV const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a*y.getRawVector() + x.getRawVector();
}

// y = a*x + y
inline void axpy(MEV& y, double const a, MEV const& x)
{
    // TODO: does that break anything?
    y.getRawVector() += a*x.getRawVector();
}

// y = a*x + y
inline void axpby(MEV& y, double const a, double const b, MEV const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a*x.getRawVector() + b*y.getRawVector();
}


// Matrix

inline void copy(MEM const& A, MEM& B)
{
    B = A;
}

// A = a*A
inline void scale(MEM& A, double const a)
{
    // TODO: does that break anything?
    A.getRawMatrix() *= a;
}

// Y = a*Y + X
inline void aypx(MEM& Y, double const a, MEM const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a*Y.getRawMatrix() + X.getRawMatrix();
}

// Y = a*X + Y
inline void axpy(MEM& Y, double const a, MEM const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a*X.getRawMatrix() + Y.getRawMatrix();
}


// Matrix and Vector

// v3 = A*v1 + v2
inline void matMult(MEM const& A, MEV const& x, MEV& y)
{
    assert(&x != &y);
    A.multiply(x, y);
}

// v3 = A*v1 + v2
inline void matMultAdd(MEM const& A, MEV const& v1, MEV const& v2, MEV& v3)
{
    assert(&v1 != &v3);
    // TODO: does that break anything?
    v3.getRawVector() = v2.getRawVector() + A.getRawMatrix()*v1.getRawVector();
}

} // namespace BLAS

} // namespace MathLib

#endif
