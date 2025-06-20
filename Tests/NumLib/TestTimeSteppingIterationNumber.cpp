/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

#include "NumLib/TimeStepping/Algorithms/IterationNumberBasedTimeStepping.h"
#include "NumLib/TimeStepping/Algorithms/MultiplyerInterpolationType.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "Tests/TestTools.h"
#include "TimeSteppingTestingTools.h"

TEST(NumLib, TimeSteppingIterationNumberBasedConstructor_emptyInput)
{
    std::vector<int> nonlinear_iteration_numbers = {};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};
    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    EXPECT_THROW(
        NumLib::IterationNumberBasedTimeStepping alg(
            1, 31, 1, 10, 1, multiplier_interpolation_type,
            std::move(nonlinear_iteration_numbers), std::move(multipliers), {}),
        std::runtime_error);
}

TEST(NumLib,
     TimeSteppingIterationNumberBasedConstructor_unsortedIterationNumbers)
{
    std::vector<int> nonlinear_iteration_numbers = {3, 2, 1};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};
    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    EXPECT_THROW(
        NumLib::IterationNumberBasedTimeStepping alg(
            1, 31, 1, 10, 1, multiplier_interpolation_type,
            std::move(nonlinear_iteration_numbers), std::move(multipliers), {}),
        std::runtime_error);
}

TEST(NumLib,
     TimeSteppingIterationNumberBasedConstructor_inputVectorSizeMismatch)
{
    std::vector<int> nonlinear_iteration_numbers = {1, 2, 3};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};
    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    EXPECT_THROW(
        NumLib::IterationNumberBasedTimeStepping alg(
            1, 31, 1, 10, 1, multiplier_interpolation_type,
            std::move(nonlinear_iteration_numbers), std::move(multipliers), {}),
        std::runtime_error);
}

TEST(NumLib, findMultiplierTimestepAcceptedPiecewiseConstant)
{
    std::vector<int> nonlinear_iteration_numbers = {1, 2, 3, 4};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};

    bool const current_time_step_is_accepted = true;
    for (auto const [n, m] :
         std::views::zip(nonlinear_iteration_numbers, multipliers))
    {
        auto const multiplier = NumLib::findMultiplier(
            n, current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers,
            NumLib::MultiplyerInterpolationType::PiecewiseConstant);

        EXPECT_EQ(m, multiplier);
    }
}

TEST(
    NumLib,
    findMultiplierTimestepAcceptedPiecewiseLinearInterpolationAtIntervalBoundaries)
{
    std::vector<int> nonlinear_iteration_numbers = {2, 4, 6, 8};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};

    bool const current_time_step_is_accepted = true;
    for (auto const [n, m] :
         std::views::zip(nonlinear_iteration_numbers, multipliers))
    {
        auto const multiplier = NumLib::findMultiplier(
            n, current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers, NumLib::MultiplyerInterpolationType::PiecewiseLinear);

        EXPECT_EQ(m, multiplier);
    }
}

TEST(
    NumLib,
    findMultiplierTimestepAcceptedPiecewiseLinearInterpolationOutsideBoundaries)
{
    std::vector<int> nonlinear_iteration_numbers = {2, 4, 6, 8};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};

    bool const current_time_step_is_accepted = true;

    {
        // only one nonlinear iteration
        auto const multiplier = NumLib::findMultiplier(
            1, current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers, NumLib::MultiplyerInterpolationType::PiecewiseLinear);
        EXPECT_EQ(multipliers.front(), multiplier);
    }

    {
        // more nonlinear iterations as in the vector
        // nonlinear_iteration_numbers
        auto const multiplier = NumLib::findMultiplier(
            nonlinear_iteration_numbers.back() + 1,
            current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers, NumLib::MultiplyerInterpolationType::PiecewiseLinear);
        EXPECT_EQ(multipliers.back(), multiplier);
    }
}

TEST(
    NumLib,
    findMultiplierTimestepAcceptedPiecewiseConstantInterpolationOutsideBoundaries)
{
    std::vector<int> nonlinear_iteration_numbers = {2, 4, 6, 8};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};

    bool const current_time_step_is_accepted = true;

    {
        // only one nonlinear iteration
        auto const multiplier = NumLib::findMultiplier(
            1, current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers,
            NumLib::MultiplyerInterpolationType::PiecewiseConstant);
        EXPECT_EQ(multipliers.front(), multiplier);
    }

    {
        // more nonlinear iterations as in the vector
        // nonlinear_iteration_numbers
        auto const multiplier = NumLib::findMultiplier(
            nonlinear_iteration_numbers.back() + 1,
            current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers, NumLib::MultiplyerInterpolationType::PiecewiseLinear);
        EXPECT_EQ(multipliers.back(), multiplier);
    }
}

TEST(NumLib,
     findMultiplierTimestepAcceptedPiecewiseLinearInterpolationInIntervalCenter)
{
    std::vector<int> nonlinear_iteration_numbers = {2, 4, 6, 8};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};

    bool const current_time_step_is_accepted = true;
    // check with iteration numbers not contained in vector
    // nonlinear_iteration_numbers
    for (std::size_t i = 1; i < nonlinear_iteration_numbers.size(); ++i)
    {
        auto const multiplier = NumLib::findMultiplier(
            (nonlinear_iteration_numbers[i] +
             nonlinear_iteration_numbers[i - 1]) /
                2,
            current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers, NumLib::MultiplyerInterpolationType::PiecewiseLinear);

        EXPECT_EQ((multipliers[i] + multipliers[i - 1]) / 2, multiplier);
    }
}

TEST(NumLib, findMultiplierTimestepRejected)
{
    std::vector<int> nonlinear_iteration_numbers = {1, 2, 3, 4};
    std::vector<double> multipliers = {2.0, 1.0, 0.5, 0.25};

    bool const current_time_step_is_accepted = false;
    for (auto const [n, m] :
         std::views::zip(nonlinear_iteration_numbers, multipliers))
    {
        auto const multiplier = NumLib::findMultiplier(
            n, current_time_step_is_accepted, nonlinear_iteration_numbers,
            multipliers,
            NumLib::MultiplyerInterpolationType::PiecewiseConstant);

        if (m < 1)
        {
            EXPECT_EQ(m, multiplier);
        }
        else
        {
            // in case of a rejected timestep and multiplier >= 1: the smallest
            // of the multipliers
            auto const min_multiplier =
                *std::min_element(multipliers.begin(), multipliers.end());

            EXPECT_EQ(min_multiplier, multiplier);
        }
    }
}

TEST(NumLib, TimeSteppingIterationNumberBased1)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    NumLib::IterationNumberBasedTimeStepping alg(
        1, 31, 1, 10, 1, multiplier_interpolation_type,
        std::move(iter_times_vector), std::move(multiplier_vector), {});

    const double solution_error = 0.;
    auto const end_time = alg.end();
    NumLib::TimeStep previous_timestep(alg.begin());
    NumLib::TimeStep current_timestep(alg.begin());
    auto [step_accepted, timestepper_dt] =
        alg.next(solution_error, 1, previous_timestep, current_timestep);
    ASSERT_TRUE(step_accepted);
    timestepper_dt = (current_timestep.current() + timestepper_dt > end_time)
                         ? end_time() - current_timestep.current()()
                         : timestepper_dt;
    NumLib::updateTimeSteps(timestepper_dt, previous_timestep,
                            current_timestep);
    ASSERT_EQ(1u, current_timestep.timeStepNumber());
    ASSERT_EQ(1., current_timestep.previous()());
    ASSERT_EQ(2., current_timestep.current()());
    ASSERT_EQ(1., current_timestep.dt());
    ASSERT_TRUE(current_timestep.isAccepted());

    auto [step_accepted1, timestepper_dt1] =
        alg.next(solution_error, 1, previous_timestep, current_timestep);
    ASSERT_TRUE(step_accepted1);
    timestepper_dt1 = (current_timestep.current() + timestepper_dt1 > end_time)
                          ? end_time() - current_timestep.current()()
                          : timestepper_dt1;
    NumLib::updateTimeSteps(timestepper_dt1, previous_timestep,
                            current_timestep);

    auto [step_accepted2, timestepper_dt2] =
        alg.next(solution_error, 3, previous_timestep, current_timestep);
    ASSERT_TRUE(step_accepted2);
    timestepper_dt2 = (current_timestep.current() + timestepper_dt2 > end_time)
                          ? end_time() - current_timestep.current()()
                          : timestepper_dt2;
    NumLib::updateTimeSteps(timestepper_dt2, previous_timestep,
                            current_timestep);
    ASSERT_EQ(3u, current_timestep.timeStepNumber());
    ASSERT_EQ(4., current_timestep.previous()());
    ASSERT_EQ(6., current_timestep.current()());
    ASSERT_EQ(2., current_timestep.dt());
    ASSERT_TRUE(current_timestep.isAccepted());

    auto [step_accepted3, timestepper_dt3] =
        alg.next(solution_error, 5, previous_timestep, current_timestep);
    ASSERT_TRUE(step_accepted3);
    timestepper_dt3 = (current_timestep.current() + timestepper_dt3 > end_time)
                          ? end_time() - current_timestep.current()()
                          : timestepper_dt3;
    NumLib::updateTimeSteps(timestepper_dt3, previous_timestep,
                            current_timestep);
    ASSERT_EQ(4u, current_timestep.timeStepNumber());
    ASSERT_EQ(6., current_timestep.previous()());
    ASSERT_EQ(7., current_timestep.current()());
    ASSERT_EQ(1., current_timestep.dt());
    ASSERT_TRUE(current_timestep.isAccepted());

    auto [step_accepted4, timestepper_dt4] =
        alg.next(solution_error, 7, previous_timestep, current_timestep);
    ASSERT_TRUE(step_accepted4);
    timestepper_dt4 = (current_timestep.current() + timestepper_dt4 > end_time)
                          ? end_time() - current_timestep.current()()
                          : timestepper_dt4;
    NumLib::updateTimeSteps(timestepper_dt4, previous_timestep,
                            current_timestep);
    ASSERT_EQ(5u, current_timestep.timeStepNumber());
    ASSERT_EQ(7., current_timestep.previous()());
    ASSERT_EQ(8., current_timestep.current()());
    ASSERT_EQ(1., current_timestep.dt());
    ASSERT_TRUE(current_timestep.isAccepted());

    auto [step_accepted5, timestepper_dt5] =
        alg.next(solution_error, 8, previous_timestep, current_timestep);
    ASSERT_TRUE(step_accepted5);
    timestepper_dt5 = (current_timestep.current() + timestepper_dt5 > end_time)
                          ? end_time() - current_timestep.current()()
                          : timestepper_dt5;
    NumLib::updateTimeSteps(timestepper_dt5, previous_timestep,
                            current_timestep);
    ASSERT_EQ(6u, current_timestep.timeStepNumber());
    ASSERT_EQ(8., current_timestep.previous()());
    ASSERT_EQ(9, current_timestep.current()());
    ASSERT_EQ(1., current_timestep.dt());
    ASSERT_TRUE(current_timestep.isAccepted());

    auto [step_accepted6, timestepper_dt6] =
        alg.next(solution_error, 4, previous_timestep, current_timestep);
    ASSERT_TRUE(step_accepted6);
    NumLib::updateTimeSteps(timestepper_dt6, previous_timestep,
                            current_timestep);
    timestepper_dt6 = (current_timestep.current() + timestepper_dt6 > end_time)
                          ? end_time() - current_timestep.current()()
                          : timestepper_dt6;
    ASSERT_EQ(7u, current_timestep.timeStepNumber());
    ASSERT_EQ(9., current_timestep.previous()());
    ASSERT_EQ(10, current_timestep.current()());
    ASSERT_EQ(1., current_timestep.dt());
    ASSERT_TRUE(current_timestep.isAccepted());
}

TEST(NumLib, TimeSteppingIterationNumberBased2)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    NumLib::IterationNumberBasedTimeStepping alg(
        1, 31, 1, 10, 1, multiplier_interpolation_type,
        std::move(iter_times_vector), std::move(multiplier_vector), {});

    std::vector<int> nr_iterations = {0, 2, 2, 2, 4, 6, 8, 4, 1};
    const std::vector<double> expected_vec_t = {1,  2,  4,  8,  16,
                                                24, 28, 29, 30, 31};

    std::vector<double> vec_t = timeStepping(alg, nr_iterations, {}, {});

    ASSERT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}

TEST(NumLib, TimeSteppingIterationNumberBased2FixedOutputTimes)
{
    std::vector<int> iter_times_vector = {0, 3, 5, 7};
    std::vector<double> multiplier_vector = {2.0, 1.0, 0.5, 0.25};
    std::vector<double> fixed_output_times = {5, 20};
    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    NumLib::IterationNumberBasedTimeStepping alg(
        1, 31, 1, 10, 1, multiplier_interpolation_type,
        std::move(iter_times_vector), std::move(multiplier_vector), {});

    std::vector<int> nr_iterations = {0, 2, 2, 2, 4, 6, 8, 4, 1, 1, 1, 1, 1};
    const std::vector<double> expected_vec_t = {1,  2,  4,  5,  7,  9,  10,
                                                11, 12, 14, 18, 20, 24, 31};

    std::vector<double> vec_t =
        timeStepping(alg, nr_iterations, fixed_output_times, {});

    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}
