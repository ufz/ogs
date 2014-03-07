/**
 * \author Norihiro Watanabe
 * \date   2012-06-25
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/Nonlinear/NewtonRaphson.h"
#include "Tests/TestTools.h"

namespace
{

typedef MathLib::DenseMatrix<double> MatrixType;
typedef MathLib::DenseVector<double> VectorType;
typedef MathLib::GaussAlgorithm<MatrixType,VectorType> DenseSolverType;

template<class F_JACOBIAN>
class ScalarDx
{
public:
    ScalarDx(F_JACOBIAN &f_J) : _f_Jacobain(f_J) {}
    // dx = - r/J
    void operator()(const double &x, const double &r, double &dx)
    {
        double j;
        _f_Jacobain(x, j);
        dx = - r/j;
    }

private:
    F_JACOBIAN &_f_Jacobain;
};

template<class F_JACOBIAN>
class VectorDx
{
public:
    VectorDx(F_JACOBIAN &f_J, MatrixType &matJ) : _f_Jacobain(f_J), _matJ(matJ) {}

    // dx = - r/J
    void operator()(const VectorType &x, const VectorType &r, VectorType &dx)
    {
        _f_Jacobain(x, _matJ);
        DenseSolverType solver(_matJ);
        VectorType rhs(r);
        rhs *= -1.;
        solver.solve(rhs, dx);
    }

private:
    F_JACOBIAN &_f_Jacobain;
    MatrixType &_matJ;
};

//##############################################################################
// Example problem 1 (one variable)
// f(x) = x*x -4 = 0
// x = 2,-2
//##############################################################################
namespace Example1
{

class Residual
{
public:
    void operator()(const double &x, double &r) { r = x*x-4.; }
};

class Jacobian
{
public:
    void operator()(const double &x, double &j) { j = 2*x; }
};

} // Example1


//##############################################################################
// Example problem 2 (two variables)
// 3x-y=-2
// 2x^2-y=0
// (x,y) = (-1/2, 1/2) and (2, 8)
//##############################################################################
namespace Example2
{

class Residual
{
public:
    void operator()(const VectorType &x, VectorType &r)
    {
        r[0] = 3*x[0]-x[1]+2.;
        r[1] = 2*x[0]*x[0]-x[1];
    }
};

class Jacobian
{
public:
    void operator()(const VectorType &x, MatrixType &j)
    {
        j(0,0) = 3.;
        j(0,1) = -1.0;
        j(1,0) = 4.*x[0];
        j(1,1) = -1.0;
    }
};

} // Example2


//##############################################################################
// Example problem 3 (10 variables)
//##############################################################################

namespace Example3
{

class Residual
{
public:
    void operator()(const VectorType &x, VectorType &r)
    {
        double P = 1.;
        double R = 10.;
        double s = sqrt(2.);
        r[1-1]= (9*P*x[1-1])/4 + (9*x[2-1]*x[3-1])/(8*s) + (P*R*x[7-1])/s;
        r[2-1]= (81*P*x[2-1])/4 + (9*x[1-1]*x[3-1])/(8*s) + (P*R*x[8-1])/s;
        r[3-1]= (-9*x[1-1]*x[2-1])/(4*s) + 9*P*x[3-1] + s*P*R*x[9-1];
        r[4-1]= 36*P*x[4-1] + s*P*R*x[10-1];
        r[5-1]= -2*x[5-1] + (x[2-1]*x[7-1])/(2*s) + (x[1-1]*x[8-1])/(2*s) - (x[4-1]*x[9-1])/s + s*x[4-1]*x[9-1] - (x[3-1]*x[10-1])/s + s*x[3-1]*x[10-1];
        r[6-1]= -8*x[6-1] - (x[1-1]*x[7-1])/s - s*x[3-1]*x[9-1];
        r[7-1]= -(x[1-1]/s) - (x[2-1]*x[5-1])/(2*s) + (x[1-1]*x[6-1])/s - (3*x[7-1])/2.0 + (3*x[3-1]*x[8-1])/(4*s) + (3*x[2-1]*x[9-1])/(4*s);
        r[8-1]= -(x[2-1]/s) - (x[1-1]*x[5-1])/(2*s) - (3*x[3-1]*x[7-1])/(4*s) - (9*x[8-1])/2.0 - (3*x[1-1]*x[9-1])/(4*s);
        r[9-1]= -(s*x[3-1]) - (x[4-1]*x[5-1])/s + s*x[3-1]*x[6-1] - (3*x[2-1]*x[7-1])/(4*s) + (3*x[1-1]*x[8-1])/(4*s) - 3*x[9-1];
        r[10-1]= -(s*x[4-1]) - (x[3-1]*x[5-1])/s - 6*x[10-1];
    }
};

class Jacobian
{
public:
    void operator()(const VectorType &x, MatrixType &j)
    {
        double P = 1.;
        double R = 10.;
        double s = sqrt(2.);
        j = .0;
        j(1-1,1-1) = (9*P)/4.0;
        j(7-1,1-1) = -(1/s)+x[6-1]/s;
        j(1-1,2-1) = (9*x[3-1])/(8*s);
        j(7-1,2-1) = -x[5-1]/(2*s) + (3*x[9-1])/(4*s);
        j(1-1,3-1) = (9*x[2-1])/(8*s);
        j(7-1,3-1) = (3*x[8-1])/(4*s);
        j(1-1,7-1) = (P*R)/s;
        j(7-1,5-1) = -x[2-1]/(2*s);
        j(2-1,1-1) = (9*x[3-1])/(8*s);
        j(7-1,6-1) = x[1-1]/s;
        j(2-1,2-1) = (81*P)/4.0;
        j(7-1,7-1) = -1.5;
        j(2-1,3-1) = (9*x[1-1])/(8*s);
        j(7-1,8-1) = (3*x[3-1])/(4*s);
        j(2-1,8-1) = (P*R)/s;
        j(7-1,9-1) = (3*x[2-1])/(4*s);
        j(3-1,1-1) = (-9*x[2-1])/(4*s);
        j(8-1,1-1) = -x[5-1]/(2*s) - (3*x[9-1])/(4*s);
        j(3-1,2-1) = (-9*x[1-1])/(4*s);
        j(8-1,2-1) = -(1/s);
        j(3-1,3-1) = 9*P;
        j(8-1,3-1) = (-3*x[7-1])/(4*s);
        j(3-1,9-1) = s*P*R;
        j(8-1,5-1) = -x[1-1]/(2*s);
        j(4-1,4-1) = 36*P;
        j(8-1,7-1) = (-3*x[3-1])/(4*s);
        j(4-1,10-1)= s*P*R;
        j(8-1,8-1) = -4.5;
        j(5-1,1-1) = x[8-1]/(2*s);
        j(8-1,9-1) = (-3*x[1-1])/(4*s);
        j(5-1,2-1) = x[7-1]/(2*s);
        j(9-1,1-1) = (3*x[8-1])/(4*s);
        j(5-1,3-1) = -(x[10-1]/s) + s*x[10-1];
        j(9-1,2-1) = (-3*x[7-1])/(4*s);
        j(5-1,4-1) = -(x[9-1]/s) + s*x[9-1];
        j(9-1,3-1) = -s + s*x[6-1];
        j(5-1,5-1) = -2.0;
        j(9-1,4-1) = -(x[5-1]/s);
        j(5-1,7-1) = x[2-1]/(2*s);
        j(9-1,5-1) = -(x[4-1]/s);
        j(5-1,8-1) = x[1-1]/(2*s);
        j(9-1,6-1) = s*x[3-1];
        j(5-1,9-1) = -(x[4-1]/s) + s*x[4-1];
        j(9-1,7-1) = (-3*x[2-1])/(4*s);
        j(5-1,10-1)= -(x[3-1]/s) + s*x[3-1];
        j(9-1,8-1) = (3*x[1-1])/(4*s);
        j(6-1,1-1) = -(x[7-1]/s);
        j(9-1,9-1) = -3.0;
        j(6-1,3-1) = -(s*x[9-1]);
        j(10-1,3-1) = -(x[5-1]/s);
        j(6-1,6-1) = -8.0;
        j(10-1,4-1) = -s;
        j(6-1,7-1) = -(x[1-1]/s);
        j(10-1,5-1) = -(x[3-1]/s);
        j(6-1,9-1) = -(s*x[3-1]);
        j(10-1,10-1)= -6.0;
    }
};

} // Example 3

} //namespace

//##############################################################################
// Tests
//##############################################################################
TEST(MathLib, NonlinearNR_double)
{
    Example1::Residual f_r;
    Example1::Jacobian f_j;
    ScalarDx<Example1::Jacobian> f_dx(f_j);
    double x0 = 6.0;
    double x = .0;
    MathLib::Nonlinear::NewtonRaphson nr;
    nr.solve(f_r, f_dx, x0, x);

    ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(MathLib, NonlinearNR_dense)
{
    Example2::Residual f_r;
    Example2::Jacobian f_j;
    MatrixType matJ(2, 2);
    VectorDx<Example2::Jacobian> f_dx(f_j, matJ);
    VectorType x0(2), x(2);
    x0 = 6.0;
    x = .0;
    MathLib::Nonlinear::NewtonRaphson nr;
    nr.solve(f_r, f_dx, x0, x);

    double my_expect[] = {2., 8.};
    ASSERT_ARRAY_NEAR(my_expect, x, 2, 1e-5);
}

TEST(MathLib, NonlinearNR_dense2)
{
    Example3::Residual f_r;
    Example3::Jacobian f_j;
    const std::size_t n = 10;
    MatrixType matJ(n, n, .0);
    VectorDx<Example3::Jacobian> f_dx(f_j, matJ);
    VectorType x0(n), x(n);
    x0 = 1.;
    x = 0.;
    MathLib::Nonlinear::NewtonRaphson nr;
    nr.solve(f_r, f_dx, x0, x);

    double my_expect[] = {3.39935, 3.70074e-018, -1.42576e-017, 1.4903e-021, 4.35602e-018, 0.325, -1.08167, -5.61495e-018, 7.58394e-018, -3.79368e-021};
    ASSERT_ARRAY_NEAR(my_expect, x, n, 1e-5);
}

