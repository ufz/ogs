// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <spdlog/logger.h>

#include <memory>
#include <string>
#include <utility>

#include "baselib_export.h"

namespace BaseLib
{
extern BASELIB_EXPORT std::shared_ptr<spdlog::logger> console;
void setConsoleLogLevel(std::string const& level_string);
void initOGSLogger(std::string const& log_level);
}  // namespace BaseLib

template <typename... Args>
void DBUG(fmt::format_string<Args...> fmt, Args&&... args)
{
    assert(BaseLib::console != nullptr);
    BaseLib::console->debug(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void INFO(fmt::format_string<Args...> fmt, Args&&... args)
{
    assert(BaseLib::console != nullptr);
    BaseLib::console->info(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void WARN(fmt::format_string<Args...> fmt, Args&&... args)
{
    assert(BaseLib::console != nullptr);
    BaseLib::console->warn(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void ERR(fmt::format_string<Args...> fmt, Args&&... args)
{
    assert(BaseLib::console != nullptr);
    BaseLib::console->error(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void CRITICAL(fmt::format_string<Args...> fmt, Args&&... args)
{
    assert(BaseLib::console != nullptr);
    BaseLib::console->critical(fmt, std::forward<Args>(args)...);
}
