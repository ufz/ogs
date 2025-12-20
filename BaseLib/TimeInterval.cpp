// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TimeInterval.h"

#include <memory>

#include "ConfigTree.h"

namespace BaseLib
{
TimeInterval createTimeInterval(ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_interval}
    auto const& time_interval_config = config.getConfigSubtree("time_interval");

    const auto start_time =
        //! \ogs_file_param{prj__time_loop__processes__process__time_interval__start}
        time_interval_config.getConfigParameter<double>("start");

    const auto end_time =
        //! \ogs_file_param{prj__time_loop__processes__process__time_interval__end}
        time_interval_config.getConfigParameter<double>("end");

    return {start_time, end_time};
}
}  // namespace BaseLib
