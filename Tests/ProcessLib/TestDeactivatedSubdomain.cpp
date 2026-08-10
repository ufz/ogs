// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ProcessLib/DeactivatedSubdomain.h"

// The time passed to the deactivation is accumulated from the time step sizes
// and may miss the interval's bounds by a few units in the last place. Such a
// time must still be inside the support interval, otherwise the deactivated
// subdomain is activated in that time step, cf. the time comparisons in the
// time loop, which are all tolerant, NumLib::Time.
TEST(ProcessLibDeactivatedSubdomain, TimeSupportIntervalToleratesRoundOff)
{
    double const start = 32105.6977705;
    double const end = 64211.395541;
    MathLib::PiecewiseLinearInterpolation const time_interval{
        std::vector<double>{start, end}, std::vector<double>{0., 1.}};

    EXPECT_TRUE(ProcessLib::isTimeInSupportInterval(time_interval, start));
    EXPECT_TRUE(ProcessLib::isTimeInSupportInterval(time_interval, end));

    // One unit in the last place outside of the interval.
    EXPECT_TRUE(ProcessLib::isTimeInSupportInterval(time_interval,
                                                    std::nextafter(start, 0.)));
    EXPECT_TRUE(ProcessLib::isTimeInSupportInterval(
        time_interval, std::nextafter(end, 2. * end)));

    // Clearly outside of the interval.
    EXPECT_FALSE(ProcessLib::isTimeInSupportInterval(time_interval, 0.));
    EXPECT_FALSE(
        ProcessLib::isTimeInSupportInterval(time_interval, start - 1.e-3));
    EXPECT_FALSE(
        ProcessLib::isTimeInSupportInterval(time_interval, end + 1.e-3));
}
