/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestNonlinearPicard.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/Nonlinear/Picard.h"
#include "Tests/TestTools.h"

namespace
{

typedef MathLib::DenseMatrix<double> MatrixType;
typedef MathLib::DenseVector<double> VectorType;
typedef MathLib::GaussAlgorithm<MatrixType,VectorType> DenseSolverType;


//##############################################################################
// Example problem 1 (one variable)
// f(x) = x*x -4 = 0
// x = 2,-2
// the above function can be transformed from Newton-Rahpson equation as follows
// x = g(x) = x - 1/f'(x) * f(x) = (x*x+4)/(2x)
//##############################################################################
class Example1
{
public:
    void operator()(const double &x, double &x_new) { x_new = (x*x+4.)/(2.*x); }
};


//##############################################################################
// Example problem 2 (two variables)
// 3x-y=-2
// 2x^2-y=0
// (x,y) = (-1/2, 1/2) and (2, 8)
//##############################################################################
class Example2
{
public:
    Example2(size_t n) : A(n, n), b(n) {}

    void operator()(VectorType &x, VectorType &x_new)
    {
        A(0,0) = 3;
        A(0,1) = -1.0;
        A(1,0) = 2*x[0];
        A(1,1) = -1;
        b[0] = -2;
        b[1] = 0.;

        DenseSolverType solver(A);
        solver.solve(b, x_new);
    }
private:
    MatrixType A;
    VectorType b;
};

} //namespace

//##############################################################################
// Tests
//##############################################################################
TEST(MathLib, NonlinearPicard_double)
{
    Example1 f1;
    double x0 = 6.0;
    double x = .0;
    MathLib::Nonlinear::Picard picard;
    //picard.printErrors(true);

    // abs tol, reach max. iterations
    picard.setMaxIterations(3);
    picard.setAbsTolerance(1e-5);
    picard.setRelTolerance(std::numeric_limits<double>::max());
    ASSERT_FALSE(picard.solve(f1, x0, x));
    ASSERT_EQ(3u, picard.getNIterations());

    // abs tol, converge
    picard.setMaxIterations(10);
    ASSERT_TRUE(picard.solve(f1, x0, x));
    ASSERT_NEAR(2.0, x, 1e-5);
    ASSERT_EQ(5u, picard.getNIterations());

    // rel tol, converge
    picard.setMaxIterations(100);
    picard.setAbsTolerance(std::numeric_limits<double>::max());
    picard.setRelTolerance(1e-5);
    ASSERT_TRUE(picard.solve(f1, x0, x));
    ASSERT_NEAR(2.0, x, 1e-5);
    ASSERT_EQ(5u, picard.getNIterations());
}

TEST(MathLib, NonlinearPicard_vector_x0)
{
    Example2 f2(2);
    VectorType x0(2), x(2);

    // initial guess1
    x0[0] = 2.;
    x0[1] = 9.;
    x = 0.0;
    MathLib::Nonlinear::Picard picard;
    ASSERT_TRUE(picard.solve(f2, x0, x));

    double my_expect1[] = {2., 8.};
    ASSERT_ARRAY_NEAR(my_expect1, x, 2, 1e-5);

    // initial guess2
    x0 = 6.0;
    x = 0.0;
    picard.solve(f2, x0, x);

    double my_expect2[] = {-0.5, 0.5};
    ASSERT_ARRAY_NEAR(my_expect2, x, 2, 1e-5);
}

TEST(MathLib, NonlinearPicard_vector_norms)
{
    Example2 f2(2);
    VectorType x0(2), x(2);
    MathLib::Nonlinear::Picard picard;
    //picard.printErrors(true);

    x0[0] = 2.;
    x0[1] = 9.;
    x = 0.0;
    double my_expect1[] = {2., 8.};

    // check relative errors with different norm types
    picard.setMaxIterations(1);
    picard.setNormType(MathLib::VecNormType::NORM1);
    ASSERT_FALSE(picard.solve(f2, x0, x));
    ASSERT_NEAR(0.1, picard.getRelError(), 1e-3);

    picard.setNormType(MathLib::VecNormType::NORM2);
    picard.setMaxIterations(1);
    ASSERT_FALSE(picard.solve(f2, x0, x));
    ASSERT_NEAR(0.1213, picard.getRelError(), 1e-3);

    picard.setNormType(MathLib::VecNormType::INFINITY_N);
    picard.setMaxIterations(1);
    ASSERT_FALSE(picard.solve(f2, x0, x));
    ASSERT_NEAR(0.125, picard.getRelError(), 1e-3);

    // solution should be converged with any norm types
    picard.setMaxIterations(5);
    picard.setNormType(MathLib::VecNormType::NORM1);
    ASSERT_TRUE(picard.solve(f2, x0, x));
    ASSERT_ARRAY_NEAR(my_expect1, x, 2, 1e-5);

    picard.setNormType(MathLib::VecNormType::NORM2);
    ASSERT_TRUE(picard.solve(f2, x0, x));
    ASSERT_ARRAY_NEAR(my_expect1, x, 2, 1e-5);

    picard.setNormType(MathLib::VecNormType::INFINITY_N);
    ASSERT_TRUE(picard.solve(f2, x0, x));
    ASSERT_ARRAY_NEAR(my_expect1, x, 2, 1e-5);

}


