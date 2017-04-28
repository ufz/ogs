/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestTimeSteppingEvolutionaryPIDcontroller.cpp
 *  Created on April 5, 2017, 12:09 PM
 */

#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"

#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TimeStepping/Algorithms/EvolutionaryPIDcontroller.h"

#include "Tests/TestTools.h"

std::unique_ptr<NumLib::ITimeStepAlgorithm> createTestTimeStepper(
    const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
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
        "   <norm_type> NORM2 </norm_type>"
        "</time_stepping>";
    auto const PIDSteper = createTestTimeStepper(xml);

    double solution_error = 0.;
    // 1st step
    ASSERT_TRUE(PIDSteper->next(solution_error));
    NumLib::TimeStep ts = PIDSteper->getTimeStep();
    double h_new = 0.01;
    double t_previous = 0.;
    ASSERT_EQ(1u, ts.steps());
    ASSERT_EQ(t_previous, ts.previous());
    ASSERT_EQ(t_previous + h_new, ts.current());
    ASSERT_EQ(h_new, ts.dt());
    ASSERT_TRUE(PIDSteper->accepted());
    t_previous += h_new;

    // e_n_minus1 is filled.
    solution_error = 1.0e-4;
    PIDSteper->next(solution_error);
    ts = PIDSteper->getTimeStep();
    h_new = ts.dt();
    ASSERT_EQ(2u, ts.steps());
    ASSERT_NEAR(t_previous, ts.previous(), 1.e-10);
    ASSERT_NEAR(t_previous + h_new, ts.current(), 1.e-10);
    ASSERT_TRUE(PIDSteper->accepted());
    t_previous += h_new;

    // e_n_minus2 is filled.
    solution_error = 0.5e-3;
    PIDSteper->next(solution_error);
    ts = PIDSteper->getTimeStep();
    h_new = ts.dt();
    ASSERT_EQ(3u, ts.steps());
    ASSERT_NEAR(t_previous, ts.previous(), 1.e-10);
    ASSERT_NEAR(t_previous + h_new, ts.current(), 1.e-10);
    ASSERT_TRUE(PIDSteper->accepted());

    // error > TOL=1.3-3, step rejected and a new step size estimated.
    solution_error = 0.01;
    PIDSteper->next(solution_error);
    ts = PIDSteper->getTimeStep();
    h_new = ts.dt();
    // No change in ts.steps
    ASSERT_EQ(3u, ts.steps());
    // No change in ts.previous(), which is the same as that of the previous
    // step.
    ASSERT_NEAR(t_previous, ts.previous(), 1.e-10);
    ASSERT_NEAR(t_previous + h_new, ts.current(), 1.e-10);
    ASSERT_FALSE(PIDSteper->accepted());
    t_previous += h_new;

    // With e_n, e_n_minus1, e_n_minus2
    solution_error = 0.4e-3;
    PIDSteper->next(solution_error);
    ts = PIDSteper->getTimeStep();
    h_new = ts.dt();
    ASSERT_EQ(4u, ts.steps());
    ASSERT_NEAR(t_previous, ts.previous(), 1.e-10);
    ASSERT_NEAR(t_previous + h_new, ts.current(), 1.e-10);
    ASSERT_TRUE(PIDSteper->accepted());
}
