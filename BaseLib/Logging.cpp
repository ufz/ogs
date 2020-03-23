/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Logging.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <map>

#include "Error.h"
namespace BaseLib
{
std::shared_ptr<spdlog::logger> console = spdlog::stdout_color_st("ogs");

void setConsoleLogLevel(std::string const& level_string)
{
    using namespace spdlog::level;
    std::map<std::string, level_enum> string_to_log_level = {
        {"none", off},  {"critical", critical}, {"error", err}, {"warn", warn},
        {"info", info}, {"debug", debug},       {"all", trace}};

    auto const level = string_to_log_level.find(level_string);
    if (level == string_to_log_level.end())
    {
        ERR("'{:s}' is not a valid log level!", level_string);
        OGS_FATAL("Wrong log level string.");
    }
    console->set_level(level->second);
    spdlog::set_default_logger(console);
}
}  // namespace BaseLib
