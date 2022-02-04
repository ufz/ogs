/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on April 5, 2017, 12:09 PM
 */

#include <gtest/gtest.h>

#include <memory>
#include <tuple>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "NumLib/TimeStepping/Algorithms/CreateEvolutionaryPIDcontroller.h"
#include "NumLib/TimeStepping/Algorithms/EvolutionaryPIDcontroller.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "Tests/TestTools.h"

std::unique_ptr<NumLib::TimeStepAlgorithm> createTestTimeStepper(
    const char xml[])
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("time_stepping");
    return NumLib::createEvolutionaryPIDcontroller(sub_config);
}

TEST(NumLibTimeStepping, testEvolutionaryPIDcontroller)
{
    const char xml[] =
        "<time_stepping>"
        "   <type>EvolutionaryPIDcontroller</type>"
        "   <t_initial> 0.0 </t_initial>"
        "   <t_end> 10 </t_end>"
        "   <dt_guess> 0.01 </dt_guess>"
        "   <dt_min> 0.001 </dt_min>"
        "   <dt_max> 1 </dt_max>"
        "   <rel_dt_min> 0.01 </rel_dt_min>"
        "   <rel_dt_max> 5 </rel_dt_max>"
        "   <tol> 1.e-3 </tol>"
        "</time_stepping>";
    auto const PIDStepper = createTestTimeStepper(xml);

    double solution_error = 0.;
    int const number_iterations = 0;
    // 1st step

    auto [step_accepted, timestepper_dt] =
        PIDStepper->next(solution_error, number_iterations);

    ASSERT_TRUE(step_accepted);
    PIDStepper->resetCurrentTimeStep(timestepper_dt);

    NumLib::TimeStep ts = PIDStepper->getTimeStep();
    double h_new = 0.01;
    double t_previous = 0.;
    ASSERT_EQ(1u, ts.timeStepNumber());
    ASSERT_EQ(t_previous, ts.previous());
    ASSERT_EQ(t_previous + h_new, ts.current());
    ASSERT_EQ(h_new, ts.dt());
    ASSERT_TRUE(PIDStepper->accepted());
    t_previous += h_new;

    // e_n_minus1 is filled.
    solution_error = 1.0e-4;
    auto [step_accepted1, timestepper_dt1] =
        PIDStepper->next(solution_error, number_iterations);
    ASSERT_TRUE(step_accepted1);
    PIDStepper->resetCurrentTimeStep(timestepper_dt1);
    ts = PIDStepper->getTimeStep();
    h_new = ts.dt();
    ASSERT_EQ(2u, ts.timeStepNumber());
    const double tol = 1.e-16;
    ASSERT_NEAR(t_previous, ts.previous(), tol);
    ASSERT_NEAR(t_previous + h_new, ts.current(), tol);
    ASSERT_TRUE(PIDStepper->accepted());
    t_previous += h_new;

    // e_n_minus2 is filled.
    solution_error = 0.5e-3;
    auto [step_accepted2, timestepper_dt2] =
        PIDStepper->next(solution_error, number_iterations);
    ASSERT_TRUE(step_accepted2);
    PIDStepper->resetCurrentTimeStep(timestepper_dt2);
    ts = PIDStepper->getTimeStep();
    h_new = ts.dt();
    ASSERT_EQ(3u, ts.timeStepNumber());
    ASSERT_NEAR(t_previous, ts.previous(), tol);
    ASSERT_NEAR(t_previous + h_new, ts.current(), tol);
    ASSERT_TRUE(PIDStepper->accepted());

    // error > TOL=1.3-3, step rejected and new step size estimated.
    solution_error = 0.01;
    auto [step_accepted3, timestepper_dt3] =
        PIDStepper->next(solution_error, number_iterations);
    ASSERT_TRUE(!step_accepted3);
    ts = PIDStepper->getTimeStep();
    h_new = ts.dt();
    // No change in ts.timeStepNumber
    ASSERT_EQ(3u, ts.timeStepNumber());
    // No change in ts.previous(), which is the same as that of the previous
    // step.
    ASSERT_NEAR(t_previous, ts.previous(), tol);
    ASSERT_NEAR(t_previous + h_new, ts.current(), tol);
    ASSERT_FALSE(PIDStepper->accepted());
    t_previous += h_new;

    // With e_n, e_n_minus1, e_n_minus2
    solution_error = 0.4e-3;
    auto [step_accepted4, timestepper_dt4] =
        PIDStepper->next(solution_error, number_iterations);
    ASSERT_TRUE(step_accepted4);
    PIDStepper->resetCurrentTimeStep(timestepper_dt4);
    ts = PIDStepper->getTimeStep();
    h_new = ts.dt();
    ASSERT_EQ(4u, ts.timeStepNumber());
    ASSERT_NEAR(t_previous, ts.previous(), tol);
    ASSERT_NEAR(t_previous + h_new, ts.current(), tol);
    ASSERT_TRUE(PIDStepper->accepted());
}
