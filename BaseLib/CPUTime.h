// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <ctime>

namespace BaseLib
{
/// Count CPU time
class CPUTime
{
public:
    /// Start the timer.
    void start() { start_time_ = clock(); }

    /// Get the elapsed time after started.
    double elapsed() const
    {
        return (clock() - start_time_) / static_cast<double>(CLOCKS_PER_SEC);
    }

private:
    double start_time_ = 0.;
};

}  // end namespace BaseLib
