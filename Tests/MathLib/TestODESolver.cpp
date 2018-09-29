/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "MathLib/ODE/ODESolverBuilder.h"

const double abs_tol = 1e-8;
const double rel_tol = 1e-8;

bool f(const double,
       MathLib::ODE::MappedConstVector<1> const& y,
       MathLib::ODE::MappedVector<1>& ydot)
{
    if (y[0] <= 0.0) return false;

    ydot[0] = -15.0 * y[0];
    return true;
}

bool df(const double /*t*/,
        MathLib::ODE::MappedConstVector<1> const& y,
        MathLib::ODE::MappedConstVector<1> const& /*ydot*/,
        MathLib::ODE::MappedMatrix<1, 1>& jac)
{
    if (y[0] <= 0.0) return false;

    jac(0, 0) = -15.0;
    return true;
}

struct ExtraData
{
    double value = 12.5;
};

bool f_extra(const double,
             MathLib::ODE::MappedConstVector<1> const& y,
             MathLib::ODE::MappedVector<1>& ydot,
             ExtraData& data)
{
    if (y[0] <= 0.0) return false;

    ydot[0] = -data.value * y[0];
    return true;
}

bool any_ode_solver_libs_available()
{
#ifdef CVODE_FOUND
    return true;
#else
    return false;
#endif  // CVODE_FOUND
}

template <unsigned NumEquations>
std::unique_ptr<MathLib::ODE::ODESolver<NumEquations>> make_ode_solver(
    boost::property_tree::ptree const& conf)
{
    // Make sure testrunner does not crash if we haven't built with support for
    // any external ODE solver lib.
    if (!any_ode_solver_libs_available())
    {
        ERR(
            "I cannot create any ODE solver. This test therefore might be "
            "skipped.");
        return nullptr;
    }

    BaseLib::ConfigTree config(conf, "", BaseLib::ConfigTree::onerror,
                               BaseLib::ConfigTree::onwarning);
    return MathLib::ODE::createODESolver<NumEquations>(config);
}

// There is no definition of this function in order to prevent passing temporary
// property trees! There will be linker errors if you do so anyway.
template <unsigned NumEquations>
std::unique_ptr<MathLib::ODE::ODESolver<NumEquations>> make_ode_solver(
    boost::property_tree::ptree&&);

void check(const double time_reached, const double y, const double y_dot,
           const double time, const double y_ana, const double y_dot_ana)
{
    DBUG("t: %14.7g, y: %14.7g, diff: %14.7g, y_dot: %14.7g, diff: %14.7g",
         time_reached, y, y - y_ana, y_dot, y_dot - y_dot_ana);
    (void)y_dot_ana;  // Avoid unused variable warning when DBUG output is
                      // disabled.

    EXPECT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());

    auto const abs_err = std::abs(y - y_ana);
    auto const rel_err = abs_err / std::abs(y_ana);
    EXPECT_GT(10.0 * abs_tol, abs_err);

    // relative errors are not checked.
    // EXPECT_GT(10.0*rel_tol, rel_err);

    // make sure that the relative error of the derivative is not much bigger
    // than the one of the solution y.
    auto const dot_rel_err = std::abs(y_dot / y_ana - 1.0);
    EXPECT_LT(10.0 * rel_err, dot_rel_err);
}

TEST(MathLibCVodeTest, Exponential)
{
    // initial values
    const double y0 = 1.0;
    const double t0 = 0.0;

    auto tree = boost::property_tree::ptree{};
    auto ode_solver = make_ode_solver<1>(tree);

    ASSERT_EQ(any_ode_solver_libs_available(), !!ode_solver);
    // Don't run the test if the ODE solver could not be constructed.
    if (!ode_solver) return;

    ode_solver->setFunction(f, nullptr);
    ode_solver->setTolerance(abs_tol, rel_tol);

    ode_solver->setIC(t0, {y0});

    ode_solver->preSolve();

    const double dt = 1e-1;

    for (unsigned i = 1; i <= 10; ++i)
    {
        const double time = dt * i;

        ASSERT_TRUE(ode_solver->solve(time));

        auto const y = ode_solver->getSolution();
        auto const time_reached = ode_solver->getTime();
        auto const y_dot = ode_solver->getYDot(time_reached, y);

        auto const y_ana = exp(-15.0 * time);
        auto const y_dot_ana = -15.0 * exp(-15.0 * time);

        check(time_reached, y[0], y_dot[0], time, y_ana, y_dot_ana);
    }
}

TEST(MathLibCVodeTest, ExponentialExtraData)
{
    // initial values
    const double y0 = 1.0;
    const double t0 = 0.0;

    auto tree = boost::property_tree::ptree{};
    auto ode_solver = make_ode_solver<1>(tree);

    ASSERT_EQ(any_ode_solver_libs_available(), !!ode_solver);
    // Don't run the test if the ODE solver could not be constructed.
    if (!ode_solver) return;

    ExtraData data;
    auto f_lambda = [&](double t,
                        MathLib::ODE::MappedConstVector<1> const& y,
                        MathLib::ODE::MappedVector<1>& ydot)
    {
        return f_extra(t, y, ydot, data);
    };

    ode_solver->setFunction(f_lambda, nullptr);

    ode_solver->setTolerance(abs_tol, rel_tol);
    ode_solver->setIC(t0, {y0});
    ode_solver->preSolve();

    const double dt = 1e-1;

    for (unsigned i = 1; i <= 10; ++i)
    {
        const double time = dt * i;

        ASSERT_TRUE(ode_solver->solve(time));

        auto const y = ode_solver->getSolution();
        auto const time_reached = ode_solver->getTime();
        auto const y_dot = ode_solver->getYDot(time_reached, y);

        auto const y_ana = exp(-data.value * time);
        auto const y_dot_ana = -data.value * exp(-data.value * time);

        check(time_reached, y[0], y_dot[0], time, y_ana, y_dot_ana);
    }

    ode_solver->setFunction(f_lambda, nullptr);
    ode_solver->preSolve();
    for (unsigned i = 11; i <= 15; ++i)
    {
        const double time = dt * i;

        ASSERT_TRUE(ode_solver->solve(time));

        auto const y = ode_solver->getSolution();
        auto const time_reached = ode_solver->getTime();
        auto const y_dot = ode_solver->getYDot(time_reached, y);

        auto const y_ana = exp(-data.value * time);
        auto const y_dot_ana = -data.value * exp(-data.value * time);

        check(time_reached, y[0], y_dot[0], time, y_ana, y_dot_ana);
    }
}

TEST(MathLibCVodeTest, ExponentialWithJacobian)
{
    // initial values
    const double y0 = 1.0;
    const double t0 = 0.0;

    auto tree = boost::property_tree::ptree{};
    auto ode_solver = make_ode_solver<1>(tree);

    ASSERT_EQ(any_ode_solver_libs_available(), !!ode_solver);
    // Don't run the test if the ODE solver could not be constructed.
    if (!ode_solver) return;

    ode_solver->setFunction(f, df);
    ode_solver->setTolerance(abs_tol, rel_tol);

    ode_solver->setIC(t0, {y0});

    ode_solver->preSolve();

    const double dt = 1e-1;

    for (unsigned i = 1; i <= 10; ++i)
    {
        const double time = dt * i;

        ASSERT_TRUE(ode_solver->solve(time));

        auto const y = ode_solver->getSolution();
        auto const time_reached = ode_solver->getTime();
        auto const y_dot = ode_solver->getYDot(time_reached, y);

        auto const y_ana = exp(-15.0 * time);
        auto const y_dot_ana = -15.0 * exp(-15.0 * time);

        check(time_reached, y[0], y_dot[0], time, y_ana, y_dot_ana);
    }
}

TEST(MathLibCVodeTest, ExponentialWithJacobianNewton)
{
    // initial values
    const double y0 = 1.0;
    const double t0 = 0.0;

    boost::property_tree::ptree tree;
    tree.put("linear_multistep_method", "BDF");
    tree.put("nonlinear_solver_iteration", "Newton");
    auto ode_solver = make_ode_solver<1>(tree);

    ASSERT_EQ(any_ode_solver_libs_available(), !!ode_solver);
    // Don't run the test if the ODE solver could not be constructed.
    if (!ode_solver) return;

    ode_solver->setFunction(f, df);
    ode_solver->setTolerance(abs_tol, rel_tol);

    ode_solver->setIC(t0, {y0});

    ode_solver->preSolve();

    const double dt = 1e-1;

    for (unsigned i = 1; i <= 10; ++i)
    {
        const double time = dt * i;

        ASSERT_TRUE(ode_solver->solve(time));

        auto const y = ode_solver->getSolution();
        auto const time_reached = ode_solver->getTime();
        auto const y_dot = ode_solver->getYDot(time_reached, y);

        auto const y_ana = exp(-15.0 * time);
        auto const y_dot_ana = -15.0 * exp(-15.0 * time);

        check(time_reached, y[0], y_dot[0], time, y_ana, y_dot_ana);
    }
}
