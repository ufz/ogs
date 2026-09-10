// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <spdlog/sinks/ostream_sink.h>

#include <memory>
#include <sstream>
#include <string>

#include "BaseLib/Logging.h"

namespace Tests
{
/// Redirects the logging facade into a string for as long as it is alive, so
/// that a warning can be asserted on. Where a correlation is capped or
/// extrapolated and the value alone no longer says so, the warning is part of
/// the behaviour under test and not a side effect of it.
class LogCapture final
{
public:
    LogCapture() : previous_(BaseLib::console)
    {
        auto sink = std::make_shared<spdlog::sinks::ostream_sink_st>(stream_);
        BaseLib::console =
            std::make_shared<spdlog::logger>("test capture", sink);
        BaseLib::console->set_level(spdlog::level::warn);
    }

    ~LogCapture() { BaseLib::console = previous_; }

    LogCapture(LogCapture const&) = delete;
    LogCapture& operator=(LogCapture const&) = delete;

    std::string text() const { return stream_.str(); }

private:
    std::shared_ptr<spdlog::logger> previous_;
    std::ostringstream stream_;
};
}  // namespace Tests
