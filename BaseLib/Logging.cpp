// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Logging.h"

#include <spdlog/common.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <exception>
#include <iostream>
#include <map>

#include "Error.h"
#include "MPI.h"

namespace
{
#ifdef USE_PETSC
static int mpi_rank = -1;
void mpi_error_handler(const std::string& msg)
{
    assert(mpi_rank != -1);
    std::cerr << "[" << mpi_rank << "] spdlog error: " << msg << std::endl;
    std::abort();
}
#endif  // USE_PETSC
void error_handler(const std::string& msg)
{
    std::cerr << "spdlog error: " << msg << std::endl;
    std::abort();
}
}  // namespace

namespace BaseLib
{
std::shared_ptr<spdlog::logger> console;

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

void initOGSLogger(std::string const& log_level)
{
    if (!console)
    {
#ifdef USE_PETSC
        console = spdlog::stdout_color_mt("ogs");
#else   // USE_PETSC
        console = spdlog::stdout_color_st("ogs");
#endif  // USE_PETSC
        // Default pattern and error handler both for MPI and non-MPI builds.
        spdlog::set_pattern("%^%l:%$ %v");
        spdlog::set_error_handler(error_handler);

#ifdef USE_PETSC
        int mpi_init;
        MPI_Initialized(&mpi_init);
        if (mpi_init == 1)
        {
            MPI_Comm_rank(BaseLib::MPI::OGS_COMM_WORLD, &mpi_rank);
            spdlog::set_pattern(fmt::format("[{}] %^%l:%$ %v", mpi_rank));
            spdlog::set_error_handler(mpi_error_handler);
        }
#endif  // USE_PETSC
    }
    BaseLib::setConsoleLogLevel(log_level);
}
}  // namespace BaseLib
