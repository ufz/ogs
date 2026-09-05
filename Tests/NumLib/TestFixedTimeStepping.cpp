// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <autocheck/autocheck.hpp>
#include <numeric>

#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"
#include "NumLib/TimeStepping/Time.h"

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

// No time steps can be generated for an empty time interval. Both
// constructors must reject it.
struct NumLibFixedTimeSteppingEmptyTimeInterval
    : public ::testing::TestWithParam<std::pair<double, double>>
{
};

TEST_P(NumLibFixedTimeSteppingEmptyTimeInterval, ConstructionThrows)
{
    auto const [t_initial, t_end] = GetParam();
    std::vector<NumLib::RepeatDtPair> const repeat_dt_pairs{{1, 1.0}};

    EXPECT_ANY_THROW(
        NumLib::FixedTimeStepping(t_initial, t_end, repeat_dt_pairs, {}));
    EXPECT_ANY_THROW(NumLib::FixedTimeStepping(t_initial, t_end, 1.0));
}

INSTANTIATE_TEST_SUITE_P(
    NumLibFixedTimeStepping, NumLibFixedTimeSteppingEmptyTimeInterval,
    // Initial time greater than the end time, initial time equal to it.
    ::testing::Values(std::pair{1.0, 0.0}, std::pair{1.0, 1.0}));

// The uniform time step size must be positive.
struct NumLibFixedTimeSteppingNonPositiveUniformDt
    : public ::testing::TestWithParam<double>
{
};

TEST_P(NumLibFixedTimeSteppingNonPositiveUniformDt, ConstructionThrows)
{
    EXPECT_ANY_THROW(NumLib::FixedTimeStepping(0.0, 1.0, GetParam()));
}

INSTANTIATE_TEST_SUITE_P(NumLibFixedTimeStepping,
                         NumLibFixedTimeSteppingNonPositiveUniformDt,
                         ::testing::Values(0.0, -1.0));

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

        NumLib::TimeStep ts_dummy(NumLib::Time(0), NumLib::Time(0), 0);
        NumLib::TimeStep ts_current(
            NumLib::Time(0), NumLib::Time(t_initial), 0);
        for (auto const& expected_time_point : expected_time_points)
        {
            auto const step_size =
                fixed_time_stepping.next(0.0 /* solution_error */,
                                         0 /* number_of_iterations */,
                                         ts_dummy,
                                         ts_current);
            // check if the current time plus the computed step size is not
            // equivalent (comparison with eps included) to the expected time
            // point
            if (ts_current.current() + step_size !=
                NumLib::Time(expected_time_point))
            {
                return false;
            }
            // check if the step size is zero
            if (step_size == 0.0)
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

    NumLib::TimeStep ts_dummy(NumLib::Time(0), NumLib::Time(0), 0);
    NumLib::TimeStep ts_current(NumLib::Time(0), NumLib::Time(t_initial), 0);
    for (auto const& expected_time_point : expected_time_points)
    {
        auto const step_size =
            fixed_time_stepping.next(0.0 /* solution_error */,
                                     0 /* number_of_iterations */,
                                     ts_dummy,
                                     ts_current);
        ASSERT_TRUE(step_size != 0.0);
        ASSERT_TRUE(ts_current.current() + step_size ==
                    NumLib::Time(expected_time_point));

        ts_current += step_size;
    }
    ASSERT_EQ(ts_current.timeStepNumber(), dts.size());
}

// If the taken time steps are smaller than the prescribed ones, the prescribed
// step sizes are used up before t_end is reached. The last prescribed step size
// must be repeated then, because a returned step size of zero stops the time
// loop.
//
// All tests of this fixture run the same time interval and take steps of the
// same size, so that they differ in the prescribed step sizes only.
struct NumLibFixedTimeSteppingMoreStepsThanPrescribed : public ::testing::Test
{
    static constexpr double t_initial = 0.0;
    static constexpr double t_end = 10.0;
    // Smaller than the prescribed step sizes of every test below, so the
    // prescribed sizes are consumed faster than the time advances. That is
    // what a process of the staggered coupling scheme experiences, since it
    // advances with the minimum step size of all processes.
    static constexpr double dt_taken = 1.0;
    // The interval is covered by steps of dt_taken, so every test lists this
    // many expected step sizes.
    static constexpr std::size_t number_of_taken_steps = 10;

    // Advances the time stepping scheme with steps of the fixed size dt_taken,
    // which is independent of the returned step sizes. The returned step sizes
    // must be the expected ones and t_end must be reached before a zero step
    // size, which stops the time loop, is returned for the first time.
    static void checkStepSizesForTakenSteps(
        NumLib::FixedTimeStepping& fixed_time_stepping,
        std::vector<double> const& expected_step_sizes)
    {
        ASSERT_EQ(number_of_taken_steps, expected_step_sizes.size());

        NumLib::TimeStep ts_dummy(NumLib::Time(0), NumLib::Time(0), 0);
        NumLib::TimeStep ts_current(NumLib::Time(0), NumLib::Time(t_initial),
                                    0);
        for (std::size_t k = 0; k < expected_step_sizes.size(); ++k)
        {
            ASSERT_EQ(expected_step_sizes[k],
                      fixed_time_stepping.next(0.0 /* solution_error */,
                                               0 /* number_of_iterations */,
                                               ts_dummy,
                                               ts_current))
                << "at time step " << k << " and time "
                << ts_current.current()();

            ts_current += dt_taken;
        }
        ASSERT_TRUE(ts_current.current() == NumLib::Time(t_end));
        ASSERT_EQ(0.0, fixed_time_stepping.next(0.0, 0, ts_dummy, ts_current));
    }
};

TEST_F(NumLibFixedTimeSteppingMoreStepsThanPrescribed, PairsNotReachingEndTime)
{
    // The prescribed steps end at 8, so the constructor appends one more step
    // of the last size 3 to reach t_end.
    std::vector<NumLib::RepeatDtPair> const repeat_dt_pairs{
        {1, 1.0}, {2, 2.0}, {1, 3.0}};
    NumLib::FixedTimeStepping fixed_time_stepping{
        t_initial, t_end, repeat_dt_pairs, {}};

    // The first five sizes are those of the scheme: the four from the pairs
    // above plus the appended one. From then on the last prescribed size of 3
    // is repeated, clamped by the remaining time in the last two steps.
    std::vector<double> const expected_step_sizes{1.0, 2.0, 2.0, 3.0, 3.0,
                                                  3.0, 3.0, 3.0, 2.0, 1.0};

    checkStepSizesForTakenSteps(fixed_time_stepping, expected_step_sizes);
}

// A fixed output time inside the last prescribed step splits that step into two
// smaller ones, so the last entry of the internal step size vector is a
// fragment of the last prescribed step size and must not be the size that is
// repeated.
TEST_F(NumLibFixedTimeSteppingMoreStepsThanPrescribed, FixedTimesForOutput)
{
    std::vector<NumLib::RepeatDtPair> const repeat_dt_pairs{{5, 2.0}};
    // 9.5 lies inside the last prescribed step [8, 10], which is therefore
    // split into 1.5 and 0.5.
    std::vector<double> const fixed_times_for_output{9.5};
    NumLib::FixedTimeStepping fixed_time_stepping{
        t_initial, t_end, repeat_dt_pairs, fixed_times_for_output};

    // The first six sizes are those of the scheme, including the split last
    // prescribed step. From then on the last prescribed size of 2 is repeated,
    // clamped by the remaining time in the final step.
    std::vector<double> const expected_step_sizes{2.0, 2.0, 2.0, 2.0, 1.5,
                                                  0.5, 2.0, 2.0, 2.0, 1.0};

    checkStepSizesForTakenSteps(fixed_time_stepping, expected_step_sizes);
}

// A (repeat, delta_t) pair whose steps all start beyond t_end is no part of the
// scheme, so its step size must not be the size that is repeated either.
TEST_F(NumLibFixedTimeSteppingMoreStepsThanPrescribed, PairBeyondEndTime)
{
    // The steps of the first two pairs reach 12 and thus already overshoot
    // t_end, so the third pair contributes no step at all.
    std::vector<NumLib::RepeatDtPair> const repeat_dt_pairs{
        {2, 1.0}, {5, 2.0}, {1, 100.0}};
    NumLib::FixedTimeStepping fixed_time_stepping{
        t_initial, t_end, repeat_dt_pairs, {}};

    // The seven sizes of the first two pairs, then their last size of 2
    // repeated---not the unused 100 of the third pair---clamped by the
    // remaining time in the final step.
    std::vector<double> const expected_step_sizes{1.0, 1.0, 2.0, 2.0, 2.0,
                                                  2.0, 2.0, 2.0, 2.0, 1.0};

    checkStepSizesForTakenSteps(fixed_time_stepping, expected_step_sizes);
}
