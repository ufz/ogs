// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
class ConfigTree;

/*!
 * Class for a time interval, which has a member to check whether the given time
 * is in this time interval.
 */
struct TimeInterval final
{
public:
    bool contains(const double current_time) const
    {
        return (current_time >= start_time && current_time <= end_time);
    }

    double start_time;
    double end_time;
};

TimeInterval createTimeInterval(ConfigTree const& config);

}  // namespace BaseLib
