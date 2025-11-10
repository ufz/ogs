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

#include <algorithm>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
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
         ranges::views::zip(nonlinear_iteration_numbers, multipliers))
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
         ranges::views::zip(nonlinear_iteration_numbers, multipliers))
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
         ranges::views::zip(nonlinear_iteration_numbers, multipliers))
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
    auto timestepper_dt =
        alg.next(solution_error, 1, previous_timestep, current_timestep);
    ASSERT_TRUE(current_timestep.isAccepted());
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

    auto timestepper_dt1 =
        alg.next(solution_error, 1, previous_timestep, current_timestep);
    ASSERT_TRUE(current_timestep.isAccepted());
    timestepper_dt1 = (current_timestep.current() + timestepper_dt1 > end_time)
                          ? end_time() - current_timestep.current()()
                          : timestepper_dt1;
    NumLib::updateTimeSteps(timestepper_dt1, previous_timestep,
                            current_timestep);

    auto timestepper_dt2 =
        alg.next(solution_error, 3, previous_timestep, current_timestep);
    ASSERT_TRUE(current_timestep.isAccepted());
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

    auto timestepper_dt3 =
        alg.next(solution_error, 5, previous_timestep, current_timestep);
    ASSERT_TRUE(current_timestep.isAccepted());
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

    auto timestepper_dt4 =
        alg.next(solution_error, 7, previous_timestep, current_timestep);
    ASSERT_TRUE(current_timestep.isAccepted());
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

    auto timestepper_dt5 =
        alg.next(solution_error, 8, previous_timestep, current_timestep);
    ASSERT_TRUE(current_timestep.isAccepted());
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

    auto timestepper_dt6 =
        alg.next(solution_error, 4, previous_timestep, current_timestep);
    ASSERT_TRUE(current_timestep.isAccepted());
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

    std::vector<int> const nr_iterations = {0, 2, 2, 2, 4, 6, 8, 4, 1};
    std::vector<double> const expected_vec_t = {1,  2,  4,  8,  16,
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

    std::vector<int> const nr_iterations = {0, 2, 2, 2, 4, 6, 8,
                                            4, 1, 1, 1, 1, 1};
    std::vector<double> const expected_vec_t = {1,  2,  4,  5,  7,  9,  10,
                                                11, 12, 14, 18, 20, 24, 31};

    std::vector<double> vec_t =
        timeStepping(alg, nr_iterations, fixed_output_times, {});

    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}

TEST(NumLib, TimeSteppingIterationNumberBased_simple)
{
    // *** initialization of IterationNumberBaseTimeStepping object
    constexpr int number_of_multipliers = 20;
    std::vector multiplier_vector(number_of_multipliers, 1.0);

    auto iter_times_vector =
        ranges::views::iota(1, static_cast<int>(multiplier_vector.size() + 1)) |
        ranges::to<std::vector>;

    std::vector<double> const fixed_output_times = {};

    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;

    double const t_initial = 0.0;
    double const t_end = 8000.0;
    double const min_dt = 0.001;
    double const max_dt = 2100;
    double const initial_dt = 100;

    NumLib::IterationNumberBasedTimeStepping alg(
        t_initial, t_end, min_dt, max_dt, initial_dt,
        multiplier_interpolation_type, std::move(iter_times_vector),
        std::move(multiplier_vector), std::move(fixed_output_times));
    // *** end initialization of IterationNumberBaseTimeStepping object

    std::vector<int> const rejected_steps = {};
    std::vector<int> const nr_iterations = {
        0, 4, 3,  9, 3,  5, 4, 3, 3, 7, 10, 4, 3, 3, 3, 5, 3, 4, 7, 15,
        3, 7, 15, 3, 11, 8, 4, 5, 5, 4, 3,  3, 3, 3, 3, 3, 3, 5, 3, 3,
        3, 3, 3,  3, 7,  5, 4, 4, 3, 3, 3,  3, 3, 3, 3, 3, 3, 3, 3, 4,
        4, 3, 3,  3, 3,  3, 3, 3, 3, 3, 3,  3, 3, 3, 3, 3, 4, 3, 3, 3};

    auto const expected_vec_t =
        ranges::views::iota(0, static_cast<int>(nr_iterations.size()) + 1) |
        ranges::views::transform([](int i) { return (i * 100.0); }) |
        ranges::to<std::vector>;

    std::vector<double> const vec_t =
        timeStepping(alg, nr_iterations, fixed_output_times, rejected_steps);

    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}

TEST(NumLib, TimeSteppingIterationNumberBased_simple2)
{
    // *** initialization of IterationNumberBaseTimeStepping object
    constexpr int number_of_multipliers = 20;
    std::vector multiplier_vector(number_of_multipliers - 1, 1.0);
    multiplier_vector.emplace_back(0.1);  // multiplier for rejected step

    auto iter_times_vector =
        ranges::views::iota(1, static_cast<int>(multiplier_vector.size() + 1)) |
        ranges::to<std::vector>;

    std::vector<double> fixed_output_times = {};

    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    double const t_initial = 0.0;
    double const t_end = 200.0;
    double const min_dt = 10;
    double const max_dt = 200;
    double const initial_dt = 100;
    NumLib::IterationNumberBasedTimeStepping alg(
        t_initial, t_end, min_dt, max_dt, initial_dt,
        multiplier_interpolation_type, std::move(iter_times_vector),
        std::move(multiplier_vector), std::move(fixed_output_times));
    // *** end initialization of IterationNumberBaseTimeStepping object

    std::vector<int> const rejected_steps = {2};
    std::vector<int> const nr_iterations = {0, 4, 3, 3,  3,  3,
                                            3, 3, 9, 10, 19, 1};
    // current time step size:
    std::vector<double> const expected_vec_t = {
        0, 100, 200, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};

    std::vector<double> const vec_t =
        timeStepping(alg, nr_iterations, fixed_output_times, rejected_steps);

    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      std::numeric_limits<double>::epsilon());
}

// extracted from Sandwich-0c-Aqeel-R5DeCoT68
TEST(NumLib, TimeSteppingIterationNumberBasedSandwich_0c)
{
    // *** initialization of IterationNumberBaseTimeStepping object
    std::vector<double> multiplier_vector = {
        5,      4.25,   3.5,    2.75,   2,      1.6333, 1.2667,
        0.9,    0.8333, 0.7667, 0.7,    0.6333, 0.5667, 0.5,
        0.4333, 0.3667, 0.3,    0.2333, 0.1667, 0.1};

    auto iter_times_vector =
        ranges::views::iota(1, static_cast<int>(multiplier_vector.size() + 1)) |
        ranges::to<std::vector>;

    std::vector<double> fixed_output_times = {
        0,       432000,  440640,  950400,  959040,  3024000, 3032640,
        5961600, 6480000, 7689600, 7776000, 7819200, 12528000};

    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    // *** original data
    double const t_initial = 0.0;
    double const t_end = 12528000.0;
    double const min_dt = 0.001;
    double const max_dt = 432000;
    double const initial_dt = 120;
    NumLib::IterationNumberBasedTimeStepping alg(
        t_initial, t_end, min_dt, max_dt, initial_dt,
        multiplier_interpolation_type, std::move(iter_times_vector),
        std::move(multiplier_vector), std::move(fixed_output_times));
    // *** end initialization of IterationNumberBaseTimeStepping object

    std::vector<int> const rejected_steps = {6, 8, 10, 12, 14, 14, 24};

    std::vector<int> const nr_iterations = {
        0,  4, 3, 9, 3, 5, 5,  // <- time step 6 rejected
        3,                     // <- time step 6 repeated
        3,  8,                 // <- time step 8 rejected
        10,                    // <- time step 8 repeated
        4,  4,                 // <- time step 10 rejected
        3,                     // <- time step 10 repeated
        3,  6,                 // <- time step 12 rejected
        3,                     // <- time step 12 repeated
        4,  8,                 // <- time step 14 rejected
        17,                    // <- time step 14 rejected again
        3,                     // <- time step 14 repeated (and accepted)
        9,  4, 5, 5, 4, 3, 3, 3, 3, 11,  // <- time step 24 rejected
        3                                // <- time step 24 repeated
    };

    // current timestep sizes:
    std::vector<double> const expected_vec_t = {
        0,
        120,
        450,
        1605,
        2567.4614999999999,
        5936.0767500000002,
        12673.30725,         // ts 6 rejected
        6609.7998000000007,  // ts 6 repeated
        8967.8304749999988,
        17220.937837499994,  // ts 8 rejected
        16395.627101249993,  // ts 8 repeated
        22090.518774595865,
        37751.470876297011,  // ts 10 rejected
        23656.61398476598,   // ts 10 repeated
        29137.947220361384,
        48322.613544945292,  // ts 12 rejected
        31056.413852819776,  // ts 12 repeated
        37771.047066424137,
        56236.288403836123,  // ts 14 rejected
        54389.764270094922,  // ts 14 rejected
        42756.662227525369,  // ts 14 repeated
        60206.315291379709,
        74747.111189489529,
        114734.29990929153,
        194708.67734889555,
        354657.43222810357,
        432000,
        440640,
        470880,
        576720,
        947160,  // ts 24 rejected
        579960,  // ts 24 repeated
        591300};

    std::vector<double> const vec_t =
        timeStepping(alg, nr_iterations, fixed_output_times, rejected_steps);

    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      5e4 * std::numeric_limits<double>::epsilon());
}

// extracted from Sandwich-0c-Aqeel-R5DeCoT68-Try47
TEST(NumLib, TimeSteppingIterationNumberBasedSandwich_0c_ChangedMultiplier)
{
    // *** initialization of IterationNumberBaseTimeStepping object
    std::vector<double> multiplier_vector = {
        2,      1.775,  1.55,   1.325,  1.1,    1.0333, 0.9667,
        0.9,    0.8667, 0.8333, 0.8,    0.7667, 0.7333, 0.7,
        0.6667, 0.6333, 0.6,    0.5667, 0.5333, 0.5};

    auto iter_times_vector =
        ranges::views::iota(1, static_cast<int>(multiplier_vector.size() + 1)) |
        ranges::to<std::vector>;

    std::vector<double> fixed_output_times = {
        0,       23759.9136,   23760,   28080,   32399.9136,
        32400,   86399.9136,   86400,   98280,   110159.9136,
        110160,  436319.9136,  436320,  440640,  444959.9136,
        444960,  522719.9136,  522720,  527040,  531359.9136,
        531360,  950399.9136,  950400,  963360,  976319.9136,
        976320,  1032479.9136, 1032480, 1036800, 1041119.9136,
        1041120, 2851199.9136, 2851200, 2864160, 2877119.9136,
        2877120, 2937599.9136, 2937600, 2948400, 2959199.9136,
        2959200, 3011039.9136, 3011040, 3026160, 3041279.9136,
        3041280, 5961599.9136, 5961600, 5963760, 5965919.9136,
        5965920, 6035039.9136, 6035040, 6041520, 6047999.9136,
        6048000, 6130079.9136, 6130080, 6132240, 6134399.9136,
        6134400, 6389279.9136, 6389280, 6391440, 6393599.9136,
        6393600, 6475679.9136, 6475680, 6477840, 6479999.9136,
        6480000, 7689599.9136, 7689600, 7691760, 7693919.9136,
        7693920, 7771679.9136, 7771680, 7788960, 7806239.9136,
        7806240, 12528000};

    NumLib::MultiplyerInterpolationType const multiplier_interpolation_type =
        NumLib::MultiplyerInterpolationType::PiecewiseConstant;
    double const t_initial = 0.0;
    double const t_end = 12528000.0;
    double const min_dt = 0.001;
    double const max_dt = 86400;
    double const initial_dt = 900;
    NumLib::IterationNumberBasedTimeStepping alg(
        t_initial, t_end, min_dt, max_dt, initial_dt,
        multiplier_interpolation_type, std::move(iter_times_vector),
        std::move(multiplier_vector), std::move(fixed_output_times));
    // *** end initialization of IterationNumberBaseTimeStepping object

    std::vector<int> const rejected_steps = {8,  8,  8,  8,  8,  14, 46, 48, 48,
                                             48, 48, 48, 48, 48, 48, 49, 49, 49,
                                             49, 49, 49, 49, 49, 49, 49, 50};

    std::vector<int> const nr_iterations = {
        0,  3,  3,  3,  3,  3,  5, 3, 8,  4, 9,  15, 4,  17, 4,  3,  3,  4,  7,
        32, 4,  3,  13, 4,  5,  4, 9, 4,  5, 4,  5,  5,  4,  4,  4,  5,  4,  3,
        3,  3,  3,  3,  3,  1,  3, 3, 3,  3, 3,  3,  3,  3,  8,  11, 4,  10, 6,
        13, 10, 10, 13, 15, 10, 1, 7, 12, 7, 19, 14, 24, 24, 15, 15, 32, 2,  36,
        1,  2,  1,  1,  1,  1,  1, 1, 1,  1, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  3,  1,  3,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  7, 3, 1,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,
        1,  3,  3,  3,  3,  3,  3, 3, 3,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  3,  3,  4,  3,  1,
        1,  1,  1,  1,  1,  1,  1, 1, 1,  1, 3,  1,  3,  3,  1,  1,  3,  3,  3,
        1,  2,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1,  2,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  1, 1, 1,  1, 1,  1,  1,  1,  1,  1,  2,  1,  1,
        3,  3,  1,  3,  3,  1,  1, 1, 1,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,
        2,  3,  3,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  4,  3,  1,  1,  1,
        1,  1,  1,  1,  1,  1,  1, 1, 3,  1, 3,  3,  3,  3,  3,  3,  3,  3,  1,
        1,  1,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  2,  2,  3,  5,  3,  3,  3,
        3,  3,  3,  1,  1,  1,  1, 1, 2,  1, 3,  2,  1,  1,  1,  1,  3,  1,  3,
        3,  1,  2,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1,  1,  2,  2,  3,  3,
        3,  3,  3,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  3, 5, 3,  3, 14, 1,  1,  1,  1,  1,  1,  1,  1,
        1,  2,  1,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  1,  1, 1, 2,  2, 2,  3,  3,  3,  3,  3,  3,  1,  1,
        1,  1,  1,  2,  1,  3,  2, 1, 1,  1, 1,  3,  3,  3,  3,  3,  1,  1,  1,
        3,  3,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1,  1,  2,  2,  2,  3,
        3,  3,  3,  3,  3,  1,  1, 1, 1,  1, 2,  2,  2,  1,  1,  1,  1,  1,  3,
        1,  3,  3,  3,  1,  1,  3, 3, 3,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  2,  2,  2, 3, 3,  3, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  3, 3, 1,  4, 3,  3,  1,  1,  1,  1,  1,  1,  1,
        3,  2,  3,  3,  3,  3,  3, 3, 3,  3, 1,  1,  1,  1,  1,  1,  1,  1,  1,
        1,  2,  2,  2,  2,  3,  3, 3, 3,  3, 3,  3,  3,  3,  3,  1,  1,  1,  1,
        2,  1,  3,  2,  1,  1,  1, 2, 3,  3, 3,  3,  3,  3,  3,  3,  1,  1,  1,
        1,  1,  1,  1,  1,  1,  1, 1, 1,  2, 2,  2,  2,  3,  3,  3,  3,  3,  3,
        3,  1,  1,  1,  2,  1,  2, 2, 2,  1, 1,  1,  2,  2,  3,  3,  1,  1,  1,
        3,  1,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1,  2,  2,  2,  2,  3,
        3,  3,  3,  3,  3,  3,  3, 3, 3,  1, 1,  1,  1,  1,  2,  1,  3,  2,  1,
        1,  2,  3,  3,  3,  1,  3, 3, 1,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,
        1,  2,  2,  2,  2,  2,  3, 3, 3,  3, 3,  3,  3,  3,  1,  1,  1,  3,  1,
        3,  2,  2,  1,  1,  1,  2, 2, 3,  3, 3,  3,  3,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  2,  2, 2, 2,  2, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  3, 3, 3,  3, 3,  4,  2,  1,  1,  1,  1,  1,  1,
        2,  1,  3,  2,  2,  3,  3, 3, 3,  3, 3,  1,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  2,  2,  2,  2,  2, 3, 3,  3, 3,  7,  3,  3,  3,  3,  3,  1,  1,
        2,  3,  3,  2,  2,  1,  1, 1, 2,  2, 2,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1,  2,  2,  2,  2,  3,
        3,  3,  7,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  3,  3,  3,  3,  3, 3, 3,  3, 3,  3,  3,  3,  3,  3,  3,  3,  3,
        3,  3,  3,  1};

    // current timestep sizes:
    std::vector<double> const expected_vec_t = {0,
                                                900,
                                                2295,
                                                4457.25,
                                                7808.7375000000002,
                                                13003.543125,
                                                21055.491843750002,
                                                23759.9136,
                                                23760,
                                                23759.99136,
                                                23759.92224,
                                                23759.921088288,
                                                23759.918592441609,
                                                23759.916096220804,
                                                23759.917593953287,
                                                23759.919578448826,
                                                23759.922654416914,
                                                23759.927422167446,
                                                23759.933739436903,
                                                23759.939846341291,
                                                23759.936792889097,
                                                23759.940838713253,
                                                23759.947109740697,
                                                23759.95170828512,
                                                23759.957801356482,
                                                23759.96450373498,
                                                23759.973384386492,
                                                23759.981081247155,
                                                23759.991279587535,
                                                23760,
                                                23760.011554546516,
                                                23760.024264547683,
                                                23760.038245548967,
                                                23760.056770375668,
                                                23760.081315771047,
                                                23760.113838419926,
                                                23760.149613333691,
                                                23760.197015094429,
                                                23760.270487823575,
                                                23760.384370553751,
                                                23760.560888785523,
                                                23760.834492044774,
                                                23761.258577096607,
                                                23761.915908926952,
                                                23763.230572587639,
                                                23765.268301261705,
                                                23768.426780706508,
                                                23773.322423845955,
                                                23780.910670712092,
                                                23792.672453354611,
                                                23810.903216450512,
                                                23839.160899249157,
                                                23882.960307587055,
                                                23878.580366753264,
                                                23910.11594075655,
                                                23951.900576310905,
                                                23944.935077563994,
                                                23927.525509160274,
                                                23922.882377267,
                                                23920.754212300708,
                                                23918.980812434296,
                                                23916.61655115784,
                                                23914.449897711089,
                                                23913.727427086767,
                                                23920.950399747202,
                                                23920.709874757609,
                                                23919.080869716003,
                                                23918.90260007645,
                                                23916.487346842165,
                                                23915.659370915546,
                                                23914.693399001157,
                                                23914.210413043962,
                                                23914.049433824428,
                                                23913.942108978765,
                                                23913.834768032768,
                                                23914.025298211916,
                                                23913.930033122342,
                                                23914.120563301491,
                                                23914.45875436948,
                                                23915.135136505454,
                                                23916.487900777411,
                                                23919.19342932132,
                                                23924.604486409135,
                                                23935.426600584768,
                                                23957.070828936037,
                                                24000.359285638573,
                                                24086.936199043641,
                                                24221.130414821499,
                                                24429.131449277182,
                                                24751.533052683488,
                                                25251.255537963261,
                                                26025.825390146907,
                                                27226.408661031564,
                                                28080,
                                                29403.066575401077,
                                                31453.81976727274,
                                                32399.9136,
                                                32400,
                                                32400.1728,
                                                32400.518400000001,
                                                32401.209600000002,
                                                32402.592000000004,
                                                32405.356800000009,
                                                32410.886400000018,
                                                32421.945600000035,
                                                32444.064000000071,
                                                32488.300800000143,
                                                32576.774400000286,
                                                32713.908480000508,
                                                32988.176640000951,
                                                33413.292288001641,
                                                34072.221542402709,
                                                35093.561886724354,
                                                36676.639420422915,
                                                39130.409597655686,
                                                42933.753372366475,
                                                48828.936223168203,
                                                57966.469641910873,
                                                72129.646440962024,
                                                86399.9136,
                                                86400,
                                                86400.083522879999,
                                                86400.212983343998,
                                                86400.471904271995,
                                                86400.989746127991,
                                                86402.025429839996,
                                                86404.096797264006,
                                                86408.239532112013,
                                                86416.525001808026,
                                                86433.095941200052,
                                                86466.237819984104,
                                                86532.521577552208,
                                                86665.089092688417,
                                                86930.224122960848,
                                                87460.494183505711,
                                                88282.412777350241,
                                                89556.386597809265,
                                                91531.046019520742,
                                                94591.768123173533,
                                                98280,
                                                103996.75940908103,
                                                110159.9136,
                                                110160,
                                                110160.1728,
                                                110160.5184,
                                                110161.2096,
                                                110162.592,
                                                110165.35680000001,
                                                110170.88640000002,
                                                110181.94560000004,
                                                110204.06400000007,
                                                110248.30080000014,
                                                110336.77440000029,
                                                110513.72160000057,
                                                110867.61600000114,
                                                111575.40480000229,
                                                112672.47744000406,
                                                114372.94003200681,
                                                117008.65704961107,
                                                121094.01842689767,
                                                127426.32856169192,
                                                137241.40927062297,
                                                152454.78436946616,
                                                176035.51577267304,
                                                212585.64944764375,
                                                269238.35664384835,
                                                355638.35664384835,
                                                436319.91359999997,
                                                436320,
                                                436320.11448000005,
                                                436320.29192400014,
                                                436320.64681200025,
                                                436321.35658800049,
                                                436322.77614000096,
                                                436325.61524400191,
                                                436331.29345200385,
                                                436342.64986800769,
                                                436365.36270001536,
                                                436410.78836403077,
                                                436501.63969206152,
                                                436683.34234812303,
                                                437046.74766024604,
                                                437610.02589403663,
                                                438736.58236161794,
                                                440482.74488636898,
                                                440640,
                                                440954.51022726204,
                                                441583.53068178613,
                                                442558.51238629845,
                                                444069.73402829259,
                                                444959.91360000003,
                                                444960,
                                                444960.15335999994,
                                                444960.46007999982,
                                                444961.07351999963,
                                                444962.30039999919,
                                                444964.75415999838,
                                                444969.66167999676,
                                                444979.47671999346,
                                                444999.10679998685,
                                                445038.3669599737,
                                                445116.8872799474,
                                                445273.92791989475,
                                                445588.00919978943,
                                                446216.17175957887,
                                                447331.16030320508,
                                                449059.39254582574,
                                                451738.15252188774,
                                                455890.23048478382,
                                                462325.95132727275,
                                                472301.31863313058,
                                                487763.13795721025,
                                                511728.95790953375,
                                                522719.91360000003,
                                                522720,
                                                522720.13391999993,
                                                522720.40175999986,
                                                522720.93743999966,
                                                522722.00879999931,
                                                522724.15151999862,
                                                522728.43695999717,
                                                522737.00783999427,
                                                522754.14959998854,
                                                522788.43311997707,
                                                522857.00015995407,
                                                522994.13423990807,
                                                523237.54723182647,
                                                523724.37321566331,
                                                524698.02518333693,
                                                526207.18573323102,
                                                527040,
                                                528705.62853353797,
                                                531287.35276052181,
                                                531359.91359999997,
                                                531360,
                                                531360.17280000006,
                                                531360.51840000018,
                                                531361.20960000041,
                                                531362.59200000088,
                                                531365.35680000181,
                                                531370.88640000368,
                                                531381.94560000743,
                                                531404.06400001491,
                                                531448.30080002989,
                                                531536.77440005983,
                                                531713.72160011972,
                                                532067.6160002395,
                                                532775.40480047907,
                                                534031.72992090427,
                                                535979.03385756339,
                                                538997.35495938489,
                                                543675.75266720844,
                                                550927.26911433483,
                                                562167.11960738071,
                                                579588.88787160197,
                                                606592.62868114468,
                                                648448.42693593609,
                                                713324.91423086275,
                                                799724.91423086275,
                                                886124.91423086275,
                                                950399.91359999997,
                                                950400,
                                                950400.11447999999,
                                                950400.29192400014,
                                                950400.6468120002,
                                                950401.35658800043,
                                                950402.7761400009,
                                                950405.61524400185,
                                                950411.29345200385,
                                                950422.64986800763,
                                                950445.36270001531,
                                                950490.78836403077,
                                                950581.63969206146,
                                                950763.34234812297,
                                                951126.74766024598,
                                                951690.02589403663,
                                                952816.58236161794,
                                                954562.74488636898,
                                                957269.29679973307,
                                                961464.45226544735,
                                                963360,
                                                966298.09898855654,
                                                970852.15242081927,
                                                976319.91359999997,
                                                976320,
                                                976320.17280000006,
                                                976320.51840000018,
                                                976321.20960000041,
                                                976322.59200000088,
                                                976325.35680000181,
                                                976330.88640000368,
                                                976341.94560000743,
                                                976364.06400001491,
                                                976408.30080002989,
                                                976496.77440005983,
                                                976673.72160011972,
                                                977027.6160002395,
                                                977735.40480047907,
                                                978991.72992090438,
                                                981221.70700965915,
                                                984678.17149722902,
                                                988480.28243355593,
                                                994373.55438486254,
                                                1003508.1259093879,
                                                1017666.7117724022,
                                                1032479.9136,
                                                1032480,
                                                1032480.1339200001,
                                                1032480.4017600002,
                                                1032480.9374400004,
                                                1032482.0088000008,
                                                1032484.1515200015,
                                                1032488.4369600029,
                                                1032496.0436160055,
                                                1032511.2569280106,
                                                1032534.8375616187,
                                                1032576.6931862728,
                                                1032660.4044355811,
                                                1032827.8269341978,
                                                1033162.671931431,
                                                1033832.3619258978,
                                                1034870.3814173212,
                                                1036800,
                                                1039790.9088031522,
                                                1041119.9136,
                                                1041120,
                                                1041120.15336,
                                                1041120.4600800001,
                                                1041121.0735200003,
                                                1041122.3004000008,
                                                1041124.7541600016,
                                                1041129.6616800033,
                                                1041139.4767200066,
                                                1041159.1068000132,
                                                1041198.3669600266,
                                                1041276.8872800531,
                                                1041433.9279201062,
                                                1041748.0092002125,
                                                1042376.1717604252,
                                                1043491.1603048026,
                                                1045470.2649710724,
                                                1048537.8772037907,
                                                1053292.6761645041,
                                                1060662.6145536096,
                                                1072086.0190567235,
                                                1089792.2960365498,
                                                1117237.0253552808,
                                                1159776.3557993136,
                                                1225712.3179875647,
                                                1312112.3179875647,
                                                1398512.3179875647,
                                                1484912.3179875647,
                                                1571312.3179875647,
                                                1657712.3179875647,
                                                1744112.3179875647,
                                                1830512.3179875647,
                                                1916912.3179875647,
                                                2003312.3179875647,
                                                2089712.3179875647,
                                                2176112.317987565,
                                                2262512.317987565,
                                                2348912.317987565,
                                                2435312.317987565,
                                                2521712.317987565,
                                                2608112.317987565,
                                                2694512.317987565,
                                                2780912.317987565,
                                                2851199.9136000001,
                                                2851200,
                                                2851200.0950399996,
                                                2851200.2423519995,
                                                2851200.4706855994,
                                                2851200.6305191191,
                                                2851200.950186159,
                                                2851201.5895202383,
                                                2851202.868188397,
                                                2851205.4255247144,
                                                2851210.5401973492,
                                                2851220.7695426187,
                                                2851241.2282331581,
                                                2851282.1456142371,
                                                2851363.9803763945,
                                                2851509.2370792236,
                                                2851799.7504848829,
                                                2852250.0462636538,
                                                2852948.0047207493,
                                                2854029.8403292475,
                                                2855706.6855224194,
                                                2858305.7955718357,
                                                2862334.4161484311,
                                                2864160,
                                                2866989.654969932,
                                                2871375.6201733262,
                                                2877119.9136000001,
                                                2877120,
                                                2877120.1727999998,
                                                2877120.5183999995,
                                                2877121.2095999988,
                                                2877122.5919999974,
                                                2877125.3567999946,
                                                2877130.886399989,
                                                2877141.9455999779,
                                                2877164.0639999555,
                                                2877208.3007999109,
                                                2877296.7743998216,
                                                2877473.7215996431,
                                                2877827.6159992861,
                                                2878535.404798572,
                                                2879791.7299173041,
                                                2882021.7070030542,
                                                2885979.9163302607,
                                                2892115.140787431,
                                                2901624.7386960443,
                                                2916364.6154543953,
                                                2937599.9136000001,
                                                2937600,
                                                2937600.1339199995,
                                                2937600.4017599993,
                                                2937600.9374399991,
                                                2937602.0087999976,
                                                2937604.1515199956,
                                                2937608.4369599912,
                                                2937616.0436159838,
                                                2937631.2569279685,
                                                2937654.8375615445,
                                                2937696.6931861425,
                                                2937780.404435338,
                                                2937947.826933729,
                                                2938282.671930511,
                                                2938952.3619240755,
                                                2939990.3814141001,
                                                2941599.3116236385,
                                                2944093.1534484229,
                                                2947958.6082768384,
                                                2948400,
                                                2949282.7834463231,
                                                2951048.3503389694,
                                                2954579.4841242619,
                                                2959199.9136000001,
                                                2959200,
                                                2959200.1727999998,
                                                2959200.5183999995,
                                                2959201.2095999988,
                                                2959202.5919999974,
                                                2959205.3567999946,
                                                2959210.886399989,
                                                2959221.9455999779,
                                                2959244.0639999555,
                                                2959288.3007999109,
                                                2959376.7743998216,
                                                2959553.7215996431,
                                                2959907.6159992861,
                                                2960615.404798572,
                                                2961871.7299173046,
                                                2964101.7070030547,
                                                2968059.9163302612,
                                                2974195.140787431,
                                                2983704.7386960443,
                                                2998444.6154543953,
                                                3011039.9136000001,
                                                3011040,
                                                3011040.1339199999,
                                                3011040.4017599998,
                                                3011040.9374399991,
                                                3011042.0087999981,
                                                3011044.1515199956,
                                                3011048.4369599917,
                                                3011056.0436159838,
                                                3011069.5454303701,
                                                3011093.5111509059,
                                                3011141.4425919778,
                                                3011237.3054741211,
                                                3011429.0312384074,
                                                3011812.4827669808,
                                                3012579.3858241267,
                                                3013768.0855627037,
                                                3016145.4850398568,
                                                3019830.4542294447,
                                                3025542.1564733055,
                                                3026160,
                                                3027395.6870533889,
                                                3029867.0611601667,
                                                3033697.6910256725,
                                                3039635.167317206,
                                                3041279.9136000001,
                                                3041280,
                                                3041280.1727999998,
                                                3041280.5183999995,
                                                3041281.2095999988,
                                                3041282.5919999974,
                                                3041285.3567999946,
                                                3041290.886399989,
                                                3041301.9455999779,
                                                3041324.0639999555,
                                                3041368.3007999109,
                                                3041456.7743998216,
                                                3041633.7215996431,
                                                3041987.6159992861,
                                                3042695.404798572,
                                                3043951.7299173046,
                                                3046181.7070030547,
                                                3050139.9163302612,
                                                3056275.140787431,
                                                3065784.7386960443,
                                                3080524.6154543953,
                                                3103371.424429839,
                                                3138783.9783417769,
                                                3193673.4369052807,
                                                3278752.0976787116,
                                                3365152.0976787116,
                                                3451552.0976787116,
                                                3537952.0976787116,
                                                3624352.0976787116,
                                                3710752.0976787116,
                                                3797152.0976787116,
                                                3883552.0976787116,
                                                3969952.0976787116,
                                                4056352.0976787116,
                                                4142752.0976787116,
                                                4229152.0976787116,
                                                4315552.0976787116,
                                                4401952.0976787116,
                                                4488352.0976787116,
                                                4574752.0976787116,
                                                4661152.0976787116,
                                                4747552.0976787116,
                                                4833952.0976787116,
                                                4920352.0976787116,
                                                5006752.0976787116,
                                                5093152.0976787116,
                                                5179552.0976787116,
                                                5265952.0976787116,
                                                5352352.0976787116,
                                                5438752.0976787116,
                                                5525152.0976787116,
                                                5611552.0976787116,
                                                5697952.0976787116,
                                                5784352.0976787116,
                                                5870752.0976787116,
                                                5957152.0976787116,
                                                5961599.9135999996,
                                                5961600,
                                                5961600.1144800009,
                                                5961600.2919240016,
                                                5961600.566962203,
                                                5961601.1170386048,
                                                5961602.2171914103,
                                                5961604.4174970193,
                                                5961608.8181082392,
                                                5961617.6193306772,
                                                5961635.2217755541,
                                                5961670.4266653089,
                                                5961724.9942444274,
                                                5961821.8516973639,
                                                5961971.9807494152,
                                                5962204.680780095,
                                                5962565.365827648,
                                                5963124.427651355,
                                                5963760,
                                                5964745.1371403998,
                                                5965919.9135999996,
                                                5965920,
                                                5965920.1728000008,
                                                5965920.5184000023,
                                                5965921.2096000053,
                                                5965922.5920000114,
                                                5965925.3568000235,
                                                5965930.8864000477,
                                                5965941.9456000961,
                                                5965964.064000193,
                                                5966008.3008003868,
                                                5966096.7744007744,
                                                5966253.8150414629,
                                                5966532.5621786835,
                                                5967027.3383472506,
                                                5967905.5660464587,
                                                5969266.8189802291,
                                                5971376.7610275745,
                                                5974647.1712009599,
                                                5979716.3069697069,
                                                5987573.4674112657,
                                                5999752.0660956809,
                                                6018628.8940565242,
                                                6035039.9135999996,
                                                6035040,
                                                6035040.1339200009,
                                                6035040.4017600017,
                                                6035040.9374400042,
                                                6035042.0088000083,
                                                6035044.1515200185,
                                                6035047.9548480352,
                                                6035055.5615040679,
                                                6035067.3518209197,
                                                6035088.2796333311,
                                                6035130.1352581549,
                                                6035213.8465078017,
                                                6035381.2690070951,
                                                6035678.4439433403,
                                                6036139.0650945222,
                                                6036853.0278788526,
                                                6037959.6701945644,
                                                6039674.9657839192,
                                                6041520,
                                                6044379.8030349249,
                                                6047999.9135999996,
                                                6048000,
                                                6048000.1728000008,
                                                6048000.5184000023,
                                                6048001.2096000053,
                                                6048002.5920000114,
                                                6048005.3568000235,
                                                6048010.8864000477,
                                                6048021.9456000961,
                                                6048044.064000193,
                                                6048088.3008003868,
                                                6048176.7744007744,
                                                6048353.7216015495,
                                                6048707.6160030998,
                                                6049335.778565852,
                                                6050450.7671147361,
                                                6052429.8717890056,
                                                6055942.7825858351,
                                                6061387.7943209196,
                                                6069827.5625103014,
                                                6082909.203203842,
                                                6103185.7462788308,
                                                6130079.9135999996,
                                                6130080,
                                                6130080.1339200009,
                                                6130080.4017600017,
                                                6130080.9374400042,
                                                6130082.0088000093,
                                                6130083.9104640177,
                                                6130087.7137920344,
                                                6130094.4646992637,
                                                6130106.447559596,
                                                6130127.7171366867,
                                                6130170.256290867,
                                                6130255.3345992276,
                                                6130425.4912159489,
                                                6130727.5192106292,
                                                6131263.6189011866,
                                                6132094.5734215518,
                                                6132240,
                                                6132530.8531568963,
                                                6133112.5594706889,
                                                6134275.9720982742,
                                                6134399.9135999996,
                                                6134400,
                                                6134400.1728000008,
                                                6134400.5184000023,
                                                6134401.2096000053,
                                                6134402.5920000114,
                                                6134405.3568000235,
                                                6134410.8864000477,
                                                6134421.9456000961,
                                                6134444.064000193,
                                                6134488.3008003868,
                                                6134576.7744007744,
                                                6134753.7216015495,
                                                6135107.6160030998,
                                                6135735.778565851,
                                                6136850.7671147361,
                                                6138829.8717890056,
                                                6142342.7825858342,
                                                6147787.7943209196,
                                                6156227.5625103004,
                                                6169309.203203842,
                                                6189585.7462788308,
                                                6221014.3880450632,
                                                6269728.7827827241,
                                                6345236.0946260979,
                                                6389279.9135999996,
                                                6389280,
                                                6389280.1339200009,
                                                6389280.4017600017,
                                                6389280.9374400042,
                                                6389282.0088000093,
                                                6389284.1515200185,
                                                6389288.4369600369,
                                                6389296.0436160704,
                                                6389311.2569281366,
                                                6389334.8375618402,
                                                6389376.693186664,
                                                6389460.4044363108,
                                                6389627.8269356033,
                                                6389925.0018718494,
                                                6390385.6230230303,
                                                6391099.5858073616,
                                                6391440,
                                                6392120.8283852767,
                                                6393176.1123824548,
                                                6393599.9135999996,
                                                6393600,
                                                6393600.1728000008,
                                                6393600.5184000023,
                                                6393601.2096000053,
                                                6393602.5920000114,
                                                6393605.3568000235,
                                                6393610.8864000477,
                                                6393621.9456000961,
                                                6393644.064000193,
                                                6393688.3008003868,
                                                6393776.7744007744,
                                                6393953.7216015495,
                                                6394267.8028829256,
                                                6394825.2971573677,
                                                6395814.8494945029,
                                                6397571.3048929172,
                                                6400689.0132251028,
                                                6405521.4611399909,
                                                6413011.7554080663,
                                                6424621.7115235841,
                                                6442617.1435026368,
                                                6470510.0630701687,
                                                6475679.9135999996,
                                                6475680,
                                                6475680.1339200009,
                                                6475680.4017600017,
                                                6475680.9374400042,
                                                6475682.0088000083,
                                                6475683.6694080159,
                                                6475686.9906240301,
                                                6475692.1385088526,
                                                6475701.2760044131,
                                                6475717.495059032,
                                                6475749.9331682706,
                                                6475814.809386746,
                                                6475944.5618236987,
                                                6476174.8723992892,
                                                6476583.6736709625,
                                                6477217.3156420561,
                                                6477840,
                                                6478805.1607548129,
                                                6479999.9135999996,
                                                6480000,
                                                6480000.1728000008,
                                                6480000.5184000023,
                                                6480001.2096000053,
                                                6480002.5920000114,
                                                6480005.3568000235,
                                                6480010.8864000477,
                                                6480021.9456000961,
                                                6480044.064000193,
                                                6480088.3008003868,
                                                6480176.7744007744,
                                                6480353.7216015495,
                                                6480667.8028829256,
                                                6481225.2971573677,
                                                6482214.8494945029,
                                                6483971.3048929172,
                                                6487089.0132251028,
                                                6491921.46113999,
                                                6499411.7554080663,
                                                6511021.7115235841,
                                                6529017.1435026368,
                                                6556910.0630701678,
                                                6600144.0883998424,
                                                6667156.8276608363,
                                                6753556.8276608363,
                                                6839956.8276608363,
                                                6926356.8276608363,
                                                7012756.8276608363,
                                                7099156.8276608363,
                                                7185556.8276608363,
                                                7271956.8276608363,
                                                7358356.8276608363,
                                                7444756.8276608363,
                                                7531156.8276608363,
                                                7617556.8276608363,
                                                7689599.9135999996,
                                                7689600,
                                                7689600.1144800009,
                                                7689600.3176820017,
                                                7689600.7240860034,
                                                7689601.5368940067,
                                                7689603.1625100141,
                                                7689606.4137420282,
                                                7689612.9162060572,
                                                7689625.9211341133,
                                                7689649.0048814146,
                                                7689695.1723760171,
                                                7689766.7319926508,
                                                7689893.7503121747,
                                                7690119.207829331,
                                                7690468.6669809222,
                                                7691010.3286658898,
                                                7691760,
                                                7692921.9905678704,
                                                7693919.9135999996,
                                                7693920,
                                                7693920.1728000008,
                                                7693920.5184000023,
                                                7693921.2096000053,
                                                7693922.5920000114,
                                                7693925.3568000235,
                                                7693930.8864000477,
                                                7693941.9456000961,
                                                7693964.064000193,
                                                7694008.3008003868,
                                                7694096.7744007744,
                                                7694253.815041462,
                                                7694532.5621786835,
                                                7695027.3383472506,
                                                7695905.5660464577,
                                                7697464.420212551,
                                                7699880.6441699946,
                                                7703625.7913040323,
                                                7709430.7693617912,
                                                7718428.4853513176,
                                                7727126.5773983933,
                                                7740608.620071359,
                                                7761505.7862144569,
                                                7771679.9135999996,
                                                7771680,
                                                7771680.1339200009,
                                                7771680.4017600017,
                                                7771680.9374400042,
                                                7771681.8882720089,
                                                7771683.3620616151,
                                                7771685.6464355048,
                                                7771689.70119916,
                                                7771696.8984046467,
                                                7771711.292815621,
                                                7771740.0816375697,
                                                7771797.659281468,
                                                7771899.8595993863,
                                                7772081.2651636908,
                                                7772403.2600403326,
                                                7772902.3520991262,
                                                7773675.9447902571,
                                                7774875.0134615097,
                                                7776733.5699019516,
                                                7779614.3323846357,
                                                7784079.5142327966,
                                                7788960,
                                                7796524.7529391656,
                                                7806239.9135999996,
                                                7806240,
                                                7806240.1728000008,
                                                7806240.5184000023,
                                                7806241.2096000053,
                                                7806242.5920000114,
                                                7806245.3568000235,
                                                7806250.8864000477,
                                                7806261.9456000961,
                                                7806284.064000193,
                                                7806328.3008003868,
                                                7806416.7744007744,
                                                7806593.7216015495,
                                                7806947.6160030998,
                                                7807575.778565852,
                                                7808690.7671147361,
                                                7810669.8717890056,
                                                7814182.7825858351,
                                                7819627.7943209196,
                                                7828067.5625103014,
                                                7841149.203203842,
                                                7853795.225262288,
                                                7873396.5594528802,
                                                7903778.6274482971,
                                                7950870.8328411933,
                                                8023863.7512001833,
                                                8110263.7512001833,
                                                8196663.7512001833,
                                                8283063.7512001833,
                                                8369463.7512001833,
                                                8455863.7512001842,
                                                8542263.7512001842,
                                                8628663.7512001842,
                                                8715063.7512001842,
                                                8801463.7512001842,
                                                8887863.7512001842,
                                                8974263.7512001842,
                                                9060663.7512001842,
                                                9147063.7512001842,
                                                9233463.7512001842,
                                                9319863.7512001842,
                                                9406263.7512001842,
                                                9492663.7512001842,
                                                9579063.7512001842,
                                                9665463.7512001842,
                                                9751863.7512001842,
                                                9838263.7512001842,
                                                9924663.7512001842,
                                                10011063.751200184,
                                                10097463.751200184,
                                                10183863.751200184,
                                                10270263.751200184,
                                                10356663.751200184,
                                                10443063.751200184,
                                                10529463.751200184,
                                                10615863.751200184,
                                                10702263.751200184,
                                                10788663.751200184,
                                                10875063.751200184,
                                                10961463.751200184,
                                                11047863.751200184,
                                                11134263.751200184,
                                                11220663.751200184,
                                                11307063.751200184,
                                                11393463.751200184,
                                                11479863.751200184,
                                                11566263.751200184,
                                                11652663.751200184,
                                                11739063.751200184,
                                                11825463.751200184,
                                                11911863.751200184,
                                                11998263.751200184,
                                                12084663.751200184,
                                                12171063.751200184,
                                                12257463.751200184,
                                                12343863.751200184,
                                                12430263.751200184,
                                                12516663.751200184,
                                                12528000.};

    std::vector<double> const vec_t =
        timeStepping(alg, nr_iterations, fixed_output_times, rejected_steps);

    ASSERT_ARRAY_NEAR(expected_vec_t, vec_t, expected_vec_t.size(),
                      5e4 * std::numeric_limits<double>::epsilon());
    EXPECT_EQ(expected_vec_t.size(), vec_t.size());
}
