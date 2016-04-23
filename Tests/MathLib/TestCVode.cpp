/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifdef CVODE_FOUND

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include "MathLib/ODE/CVodeSolver.h"
#include "MathLib/ODE/OdeSolver.h"
#include "MathLib/ODE/ConcreteOdeSolver.h"

bool f(const double,
       MathLib::MappedConstVector<1> const y,
       MathLib::MappedVector<1> ydot)
{
	if (y[0] <= 0.0) return false;

	ydot[0] = -15.0 * y[0];
	return true;
}

bool df(const double /*t*/,
        MathLib::MappedConstVector<1> const y,
        MathLib::MappedConstVector<1> /*ydot*/,
        MathLib::MappedMatrix<1, 1> jac)
{
	if (y[0] <= 0.0) return false;

	jac(0, 0) = -15.0;
	return true;
}

struct ExtraData
{
	double value = 15.0;
};

bool f_extra(const double,
             MathLib::MappedConstVector<1> const y,
             MathLib::MappedVector<1> ydot,
             ExtraData& data)
{
	if (y[0] <= 0.0) return false;

	ydot[0] = -data.value * y[0];
	return true;
}

TEST(MathLibCVodeTest, Exponential)
{
	// initial values
	const double y0 = 1.0;
	const double t0 = 0.0;

	auto tree = boost::property_tree::ptree{};
	BaseLib::ConfigTree config(tree, "",
	                           BaseLib::ConfigTree::onerror,
	                           BaseLib::ConfigTree::onwarning);
	auto ode_solver = MathLib::createOdeSolver<1>(config);

	ode_solver->setFunction(f, nullptr);
	ode_solver->setTolerance(1e-8, 1e-6);

	ode_solver->setIC(t0, {y0});

	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ASSERT_TRUE(ode_solver->solve(time));

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();
		auto const y_dot = ode_solver->getYDot(time_reached, y);

		DBUG(
		    "t: %14.7g, y: %14.7g, diff: %14.7g, y_dot: %14.7g, diff: %14.7g",
		    time_reached, y[0], y[0] - exp(-15.0 * time_reached), y_dot[0],
		    y_dot[0] + 15.0 * exp(-15.0 * time_reached));

		EXPECT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
	}
}

TEST(MathLibCVodeTest, ExponentialExtraData)
{
	// initial values
	const double y0 = 1.0;
	const double t0 = 0.0;

	auto tree = boost::property_tree::ptree{};
	BaseLib::ConfigTree config(tree, "",
	                           BaseLib::ConfigTree::onerror,
	                           BaseLib::ConfigTree::onwarning);
	auto ode_solver = MathLib::createOdeSolver<1>(config);

	ExtraData data;
	auto f_lambda = [&](double t,
		MathLib::MappedConstVector<1> const y,
		MathLib::MappedVector<1> ydot)
	{
		return f_extra(t, y, ydot, data);
	};

	ode_solver->setFunction(f_lambda, nullptr);

	ode_solver->setTolerance(1e-8, 1e-6);
	ode_solver->setIC(t0, {y0});
	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ASSERT_TRUE(ode_solver->solve(time));

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();
		auto const y_dot = ode_solver->getYDot(time_reached, y);

		DBUG(
		    "t: %14.7g, y: %14.7g, diff: %14.7g, y_dot: %14.7g, diff: %14.7g",
		    time_reached, y[0], y[0] - exp(-15.0 * time_reached), y_dot[0],
		    y_dot[0] + 15.0 * exp(-15.0 * time_reached));

		EXPECT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
	}

	ode_solver->setFunction(f_lambda, nullptr);
	ode_solver->preSolve();
	for (unsigned i = 11; i <= 15; ++i)
	{
		const double time = dt * i;

		ASSERT_TRUE(ode_solver->solve(time));

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();
		auto const y_dot = ode_solver->getYDot(time_reached, y);

		DBUG(
		    "t: %14.7g, y: %14.7g, diff: %14.7g, y_dot: %14.7g, diff: %14.7g",
		    time_reached, y[0], y[0] - exp(-15.0 * time_reached), y_dot[0],
		    y_dot[0] + 15.0 * exp(-15.0 * time_reached));

		EXPECT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
	}
}

TEST(MathLibCVodeTest, ExponentialWithJacobian)
{
	// initial values
	const double y0 = 1.0;
	const double t0 = 0.0;

	auto tree = boost::property_tree::ptree{};
	BaseLib::ConfigTree config(tree, "",
	                           BaseLib::ConfigTree::onerror,
	                           BaseLib::ConfigTree::onwarning);
	auto ode_solver = MathLib::createOdeSolver<1>(config);

	ode_solver->setFunction(f, df);
	ode_solver->setTolerance(1e-10, 1e-8);

	ode_solver->setIC(t0, {y0});

	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ASSERT_TRUE(ode_solver->solve(time));

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();
		auto const y_dot = ode_solver->getYDot(time_reached, y);

		DBUG(
		    "t: %14.7g, y: %14.7g, diff: %14.7g, y_dot: %14.7g, diff: %14.7g",
		    time_reached, y[0], y[0] - exp(-15.0 * time_reached), y_dot[0],
		    y_dot[0] + 15.0 * exp(-15.0 * time_reached));

		EXPECT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
	}
}

TEST(MathLibCVodeTest, ExponentialWithJacobianNewton)
{
	// initial values
	const double y0 = 1.0;
	const double t0 = 0.0;

	boost::property_tree::ptree conf;
	conf.put("linear_multistep_method", "BDF");
	conf.put("nonlinear_solver_iteration", "Newton");
	BaseLib::ConfigTree config(conf, "", BaseLib::ConfigTree::onerror,
	                           BaseLib::ConfigTree::onwarning);

	auto ode_solver = MathLib::createOdeSolver<1>(config);

	ode_solver->setFunction(f, df);
	ode_solver->setTolerance(1e-6, 1e-6);

	ode_solver->setIC(t0, {y0});

	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ASSERT_TRUE(ode_solver->solve(time));

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();
		auto const y_dot = ode_solver->getYDot(time_reached, y);

		DBUG(
		    "t: %14.7g, y: %14.7g, diff: %14.7g, y_dot: %14.7g, diff: %14.7g",
		    time_reached, y[0], y[0] - exp(-15.0 * time_reached), y_dot[0],
		    y_dot[0] + 15.0 * exp(-15.0 * time_reached));

		EXPECT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
	}
}

#endif  // CVODE_FOUND
