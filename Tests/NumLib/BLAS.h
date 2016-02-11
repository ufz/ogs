#pragma once

#include<cassert>
#include<Eigen/SparseCore>

namespace BLAS
{

/*
template<class Mat>
void matScale(Mat& A, double const a);
*/
















using ESM = Eigen::SparseMatrix<double>;
using EV = Eigen::VectorXd;


// Vector

void copy(EV const& x, EV& y)
{
    y = x;
}

void scale(EV& x, double const a)
{
    x *= a;
}

// y = a*y + X
void aypx(EV& y, double const a, EV const& x)
{
    y = a*y + x;
}

// y = a*x + y
void axpy(EV& y, double const a, EV const& x)
{
    y += a*x;
}


// Matrix

void copy(ESM const& A, ESM& B)
{
    B = A;
}

// A = a*A
void scale(ESM& A, double const a)
{
    A *= a;
}

// Y = a*Y + X
void aypx(ESM& Y, double const a, ESM const& X)
{
    Y = a*Y + X;
}

// Y = a*X + Y
void axpy(ESM& Y, double const a, ESM const& X)
{
    Y = a*X + Y;
}


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(ESM const& A, EV const& x, EV& y)
{
    assert(&x != &y);
    y = A*x;
}

// v3 = A*v1 + v2
void matMultAdd(ESM const& A, EV const& v1, EV const& v2, EV& v3)
{
    assert(&v1 != &v3);
    v3 = v2 + A*v1;
}


} // namespace BLAS
