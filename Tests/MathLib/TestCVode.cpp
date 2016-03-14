/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <array>

#ifdef ZEOLITE
#include "MaterialsLib/adsorption/adsorption.h"
#endif

#include "MathLib/ODE/CVodeSolver.h"
#include "MathLib/ODE/OdeSolver.h"
#include "MathLib/ODE/OdeSolverFactory.h"

#include <cstdio>

#ifdef ZEOLITE
using namespace Ads;

const unsigned NEQ = 2;  // number of equations

const double T = 303.15;  // K
const double R = 8.314;   // J/mol/K
const double M = 0.018;   // kg/mol
const double phi = 0.4;
const double rho_SR0 = 1160.0;  // kg/m^3
// const double k = 6e-3;   // s^-1

Adsorption* ads;

bool f_zeolite(const double, double const* const y, double* const ydot)
{
	const double pV = y[0];
	const double C = y[1];

	// std::printf("pV: %14.7g, C: %14.7g\n", pV, C);

	// if (0.0 < pV && pV < 1e-6) pV = 1e-6;

	if (pV >= 1.0)
	{
		ydot[1] = ads->get_reaction_rate(pV, T, M, C);
		ydot[0] = -R * T / M / phi * (1.0 - phi) * rho_SR0 * ydot[1];
		return true;
	}
	else if (pV > 0.0)
	{
		ydot[0] = ydot[1] = 0.0;
		return true;
	}
	else
	{
		ydot[0] = ydot[1] = 0.0;
		return false;
	}
}

#endif

bool f(const double,
       BaseLib::ArrayRef<const double, 1> y,
       BaseLib::ArrayRef<double, 1> ydot)
{
	if (y[0] <= 0.0) return false;

	ydot[0] = -15.0 * y[0];
	return true;
}

bool df(const double /*t*/,
        BaseLib::ArrayRef<const double, 1> y,
        BaseLib::ArrayRef<const double, 1> /*ydot*/,
        BaseLib::MatrixRef<double, 1, 1> jac)
{
	if (y[0] <= 0.0) return false;

	jac(0, 0) = -15.0;

	// std::printf("jac: %g\n", *jac);
	return true;
}

struct ExtraData
{
	double value = 15.0;
};

bool f_extra(const double,
             BaseLib::ArrayRef<const double, 1> y,
             BaseLib::ArrayRef<double, 1> ydot,
             ExtraData& data)
{
	if (y[0] <= 0.0) return false;

	ydot[0] = -data.value * y[0];
	return true;
}

#ifdef ZEOLITE
TEST(MathLibCVodeTest, ZeoliteAdsorption)
{
	// initial values
	const double C = 0.0;
	const double pV = 100.0;  // Pa

	ads = Adsorption::newInstance(SolidReactiveSystem::Z13XBF_Hauer);

	MathLib::CVodeSolverInternal::ConfigTree config;
	config.put("linear_multistep_method", "BDF");
	auto ode_solver = MathLib::createOdeSolver<2>(config);

	ode_solver->init();
	ode_solver->setTolerance({1.0, 0.01}, 1.0);

	ode_solver->setFunction(f_zeolite, nullptr);

	ode_solver->setIC(0.0, {pV, C});

	ode_solver->preSolve();

	const double dt = 1e-4;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ode_solver->solve(time);

		double const* y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();

		// std::printf("t: %g, pV: %g, C: %g\n", time_reached, y[0], y[1]);
		std::printf("[ %23.16g, %23.16g, %13.16g ],\n", time_reached, y[0],
		            y[1]);

		ASSERT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
		// std::printf("time: %g\n", time_reached);
	}
}
#endif

TEST(MathLibCVodeTest, Exponential)
{
	// initial values
	const double y0 = 1.0;
	const double t0 = 0.0;

	BaseLib::ConfigTree config(boost::property_tree::ptree{}, "",
	                           BaseLib::ConfigTree::onerror,
	                           BaseLib::ConfigTree::onwarning);
	auto ode_solver = MathLib::createOdeSolver<1>(config);

	ode_solver->init();
	ode_solver->setTolerance(1e-8, 1e-6);

	ode_solver->setFunction(f, nullptr);

	ode_solver->setIC(t0, {y0});

	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ode_solver->solve(time);

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();

		std::printf("t: %14.7g, y: %14.7g, diff: %14.7g\n", time_reached, y[0],
		            y[0] - exp(-15.0 * time_reached));

		ASSERT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
		// std::printf("time: %g\n", time_reached);
	}
}

TEST(MathLibCVodeTest, ExponentialExtraData)
{
	// initial values
	const double y0 = 1.0;
	const double t0 = 0.0;

	BaseLib::ConfigTree config(boost::property_tree::ptree{}, "",
	                           BaseLib::ConfigTree::onerror,
	                           BaseLib::ConfigTree::onwarning);
	auto ode_solver = MathLib::createOdeSolver<1, ExtraData>(config);

	ode_solver->init();
	ode_solver->setTolerance(1e-8, 1e-6);

	ExtraData data;
	ode_solver->setFunction(f_extra, nullptr, &data);

	ode_solver->setIC(t0, {y0});

	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ode_solver->solve(time);

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();

		std::printf("t: %14.7g, y: %14.7g, diff: %14.7g\n", time_reached, y[0],
		            y[0] - exp(-15.0 * time_reached));

		ASSERT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
		// std::printf("time: %g\n", time_reached);
	}
}

TEST(MathLibCVodeTest, ExponentialWithJacobian)
{
	// initial values
	const double y0 = 1.0;
	const double t0 = 0.0;

	BaseLib::ConfigTree config(boost::property_tree::ptree{}, "",
	                           BaseLib::ConfigTree::onerror,
	                           BaseLib::ConfigTree::onwarning);
	auto ode_solver = MathLib::createOdeSolver<1>(config);

	ode_solver->init();
	ode_solver->setTolerance(1e-10, 1e-8);

	ode_solver->setFunction(f, df);

	ode_solver->setIC(t0, {y0});

	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ode_solver->solve(time);

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();

		std::printf("t: %14.7g, y: %14.7g, diff: %14.7g\n", time_reached, y[0],
		            y[0] - exp(-15.0 * time_reached));

		ASSERT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
		// std::printf("time: %g\n", time_reached);
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

	ode_solver->init();
	ode_solver->setTolerance(1e-6, 1e-6);

	ode_solver->setFunction(f, df);

	ode_solver->setIC(t0, {y0});

	ode_solver->preSolve();

	const double dt = 1e-1;

	for (unsigned i = 1; i <= 10; ++i)
	{
		const double time = dt * i;

		ode_solver->solve(time);

		auto const y = ode_solver->getSolution();
		double time_reached = ode_solver->getTime();

		std::printf("t: %14.7g, y: %14.7g, diff: %14.7g\n", time_reached, y[0],
		            y[0] - exp(-15.0 * time_reached));

		ASSERT_NEAR(time, time_reached, std::numeric_limits<double>::epsilon());
		// std::printf("time: %g\n", time_reached);
	}
}
