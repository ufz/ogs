/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    BaseLib::console->debug(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void INFO(fmt::format_string<Args...> fmt, Args&&... args)
{
    BaseLib::console->info(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void WARN(fmt::format_string<Args...> fmt, Args&&... args)
{
    BaseLib::console->warn(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void ERR(fmt::format_string<Args...> fmt, Args&&... args)
{
    BaseLib::console->error(fmt, std::forward<Args>(args)...);
}
template <typename... Args>
void CRITICAL(fmt::format_string<Args...> fmt, Args&&... args)
{
    BaseLib::console->critical(fmt, std::forward<Args>(args)...);
}
