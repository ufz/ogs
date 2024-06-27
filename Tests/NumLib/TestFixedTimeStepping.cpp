/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <autocheck/autocheck.hpp>
#include <numeric>

#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

namespace ac = autocheck;

class NumLibFixedTimeStepping : public ::testing::Test
{
public:
    NumLibFixedTimeStepping()
    {
        double_classifier.trivial([](const std::vector<double>& time_step_sizes)
                                  { return time_step_sizes.empty(); });
        double_classifier.collect(
            [](const std::vector<double>& time_step_sizes)
            {
                return time_step_sizes.size() == 1
                           ? "of test cases with one time step size"
                           : "of test cases with more than one time step size "
                             "entry";
            });

        pair_classifier.trivial(
            [](const std::vector<NumLib::RepeatDtPair>& pairs)
            { return pairs.size() == 1; });
        pair_classifier.collect(
            [](const std::vector<NumLib::RepeatDtPair>& pairs)
            {
                return pairs.size() == 1
                           ? "of test cases with one pair of RepeatDtPair"
                           : "of test cases with more than one pair of "
                             "RepeatDtPair";
            });
    }

protected:
    ac::gtest_reporter gtest_reporter;
    ac::classifier<std::vector<double>> double_classifier;
    ac::classifier<std::vector<NumLib::RepeatDtPair>> pair_classifier;
};

TEST_F(NumLibFixedTimeStepping, EmptyRepeatDtPairs)
{
    std::vector<NumLib::RepeatDtPair> empty;
    ASSERT_FALSE(NumLib::FixedTimeStepping::areRepeatDtPairsValid(empty));
}

TEST_F(NumLibFixedTimeStepping, RepeatZeroDtPairs)
{
    std::vector<NumLib::RepeatDtPair> zero_repeat_dt_pair_vec{{0, 1.0}};
    ASSERT_FALSE(NumLib::FixedTimeStepping::areRepeatDtPairsValid(
        zero_repeat_dt_pair_vec));
}

std::vector<double> transformTimesToDts(std::vector<double> const& times)
{
    std::vector<double> dts(times.size());
    std::adjacent_difference(times.begin(), times.end(), dts.begin());
    return dts;
}

std::vector<NumLib::RepeatDtPair> transformToRepeatDtPair(
    std::vector<double> const& dts)
{
    std::vector<NumLib::RepeatDtPair> repeat_dt_pairs;
    std::transform(dts.begin(),
                   dts.end(),
                   std::back_inserter(repeat_dt_pairs),
                   [](auto const dt) { return std::tuple(1, dt); });
    return repeat_dt_pairs;
}

TEST_F(NumLibFixedTimeStepping, next)
{
    auto test = [](std::vector<double>& expected_time_points) -> bool
    {
        double const t_initial = expected_time_points.front();

        expected_time_points.erase(expected_time_points.begin());

        double const t_end = expected_time_points.back();
        auto dts = transformTimesToDts(expected_time_points);
        dts.front() -= t_initial;
        auto const repeat_dt_pair = transformToRepeatDtPair(dts);
        NumLib::FixedTimeStepping fixed_time_stepping{
            t_initial, t_end, repeat_dt_pair, {}};

        NumLib::TimeStep ts_dummy(0, 0, 0);
        NumLib::TimeStep ts_current(0, t_initial, 0);
        for (auto const& expected_time_point : expected_time_points)
        {
            auto [is_next, step_size] =
                fixed_time_stepping.next(0.0 /* solution_error */,
                                         0 /* number_of_iterations */,
                                         ts_dummy,
                                         ts_current);
            // this only happens if the last time step was processed or the
            // current time is already at the end time up to machine precision
            if (!is_next && step_size != 0.0)
            {
                return false;
            }
            // if the current time plus the computed step size minus the
            // expected time is larger than the minimal time step size then the
            // step size should be larger
            // if next is true then the step
            if (is_next &&
                (ts_current.current() + step_size) != expected_time_point)
            {
                return false;
            }
            // if next is true then the step size should be larger than zero
            if (is_next && step_size == 0.0)
            {
                return false;
            }
            ts_current += step_size;
        }
        return ts_current.timeStepNumber() == dts.size();
    };

    auto ordered_list_generator = ac::ordered_list(ac::generator<double>());
    auto time_points = ac::make_arbitrary(ordered_list_generator);
    // generated list must not be empty
    time_points.discard_if([](std::vector<double> const& xs)
                           { return xs.size() <= 1; });

    ac::check<std::vector<double>>(
        test, 1000, time_points, gtest_reporter, double_classifier);
}

TEST_F(NumLibFixedTimeStepping, next_StaticTest)
{
    std::vector<double> expected_time_points{};
    for (int i = 0; i < 101; ++i)
    {
        expected_time_points.push_back(i * 1e-2);
    }

    std::vector<double> fixed_output_times{};
    for (int i = 0; i < 10; ++i)
    {
        fixed_output_times.push_back(i * 1e-1);
    }

    double const t_initial = expected_time_points.front();

    expected_time_points.erase(expected_time_points.begin());

    double const t_end = expected_time_points.back();
    auto dts = transformTimesToDts(expected_time_points);
    dts.front() -= t_initial;
    auto const repeat_dt_pair = transformToRepeatDtPair(dts);
    NumLib::FixedTimeStepping fixed_time_stepping{
        t_initial, t_end, repeat_dt_pair, fixed_output_times};

    NumLib::TimeStep ts_dummy(0, 0, 0);
    NumLib::TimeStep ts_current(0, t_initial, 0);
    for (auto const& expected_time_point : expected_time_points)
    {
        auto [is_next, step_size] =
            fixed_time_stepping.next(0.0 /* solution_error */,
                                     0 /* number_of_iterations */,
                                     ts_dummy,
                                     ts_current);
        // this only happens if the last time step was processed or the
        // current time is already at the end time up to machine precision
        ASSERT_FALSE(!is_next && step_size != 0.0);

        // if the current time plus the computed step size minus the
        // expected time is larger than the minimal time step size then the
        // step size should be larger
        // if next is true then the step
        ASSERT_FALSE(is_next &&
                     (ts_current.current() + step_size) != expected_time_point);

        // if next is true then the step size should be larger than zero
        ASSERT_FALSE(is_next && step_size == 0.0);
        ts_current += step_size;
    }
    ASSERT_EQ(ts_current.timeStepNumber(), dts.size());
}
