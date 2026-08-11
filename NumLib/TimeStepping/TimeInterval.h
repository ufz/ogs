// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
/*!
 * Class for a time interval, which has a member to check whether the given time
 * is in this time interval.
 */
struct TimeInterval final
{
    /// \returns true if the given time is inside the closed interval.
    ///
    /// The comparison is done with the same tolerance as all other time
    /// comparisons in the time loop, cf. NumLib::Time, because the current time
    /// is an accumulated sum of time step sizes and may miss the interval's end
    /// by a few units in the last place.
    bool contains(double const current_time) const;

    double start_time;
    double end_time;
};

TimeInterval createTimeInterval(BaseLib::ConfigTree const& config);

}  // namespace NumLib
