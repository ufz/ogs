/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   TimeInterval.h
 *
 * Created on November 26, 2018, 4:44 PM
 */
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
