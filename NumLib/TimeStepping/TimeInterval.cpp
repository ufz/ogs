/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   TimeInterval.cpp
 *
 * Created on November 27, 2018, 5:06 PM
 *
 */

#include "TimeInterval.h"

#include "BaseLib/ConfigTree.h"

namespace NumLib
{
std::unique_ptr<TimeInterval> createTimeInterval(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_interval}
    auto const& time_interval_config = config.getConfigSubtree("time_interval");

    const double start_time =
        //! \ogs_file_param{prj__time_loop__processes__process__time_interval__start_time}
        time_interval_config.getConfigParameter<double>("start_time");

    const double end_time =
        //! \ogs_file_param{prj__time_loop__processes__process__time_interval__end_time}
        time_interval_config.getConfigParameter<double>("end_time");

    return std::make_unique<NumLib::TimeInterval>(start_time, end_time);
}
}  // end of namespace
