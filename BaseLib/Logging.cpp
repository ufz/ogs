// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Logging.h"

#include <spdlog/common.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#if defined(_WIN32) && defined(OGS_BUILD_WHEEL)
#include <io.h>
#include <spdlog/details/null_mutex.h>
#include <spdlog/sinks/base_sink.h>

#include <cstdio>
#include <mutex>
#endif

#include <exception>
#include <iostream>
#include <map>

#include "Error.h"
#include "MPI.h"

namespace
{
#if defined(_WIN32) && defined(OGS_BUILD_WHEEL)
// Stdout sink for Windows wheel builds resolving the target on every write.
//
// spdlog's built-in Windows stdout sinks (wincolor_sink and
// stdout_sink_base) cache the stdout OS HANDLE once at sink creation and
// WriteFile() to that value forever. When OGS runs in-process (e.g. via the
// pybind wheel) and the host redirects fd 1/2 afterwards (pytest fd capture
// does dup2 per test), the CRT closes the previously cached handle. Windows
// recycles handle values, so subsequent console writes can land in an
// arbitrary file OGS opened later. Standalone executables are unaffected -
// their std handles are fixed at process startup - so the stock colour sink
// is kept there.
//
// _write(_fileno(stdout), ...) goes through the live CRT fd table at call
// time, so it always follows the current redirection target. spdlog's
// formatter already appends the platform eol, which survives the raw
// _write. Trade-off: no Windows console colour.
//
// Another option would be running the wheel pytests in a subprocess but that
// hides the defect.
template <typename Mutex>
class Fd1StdoutSink final : public spdlog::sinks::base_sink<Mutex>
{
protected:
    void sink_it_(spdlog::details::log_msg const& msg) override
    {
        spdlog::memory_buf_t formatted;
        this->formatter_->format(msg, formatted);
        _write(_fileno(stdout), formatted.data(),
               static_cast<unsigned int>(formatted.size()));
    }

    void flush_() override
    {
        // _write is unbuffered; nothing to flush.
    }
};
#endif  // defined(_WIN32) && defined(OGS_BUILD_WHEEL)

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
#if defined(_WIN32) && defined(OGS_BUILD_WHEEL)
        // Custom sink re-resolving stdout per write; see Fd1StdoutSink above.
        // No colour on Windows console as a consequence.
#if defined(USE_PETSC) || defined(_OPENMP)
        console = std::make_shared<spdlog::logger>(
            "ogs", std::make_shared<Fd1StdoutSink<std::mutex>>());
#else
        console = std::make_shared<spdlog::logger>(
            "ogs",
            std::make_shared<Fd1StdoutSink<spdlog::details::null_mutex>>());
#endif
        spdlog::register_logger(console);
#else
#if defined(USE_PETSC) || defined(_OPENMP)
        // Thread-safe sink is required when OGS logs concurrently, e.g. from
        // OpenMP assembly threads. Otherwise the non-thread-safe sink races and
        // may corrupt output (e.g. console output interleaved into other
        // files).
        console = spdlog::stdout_color_mt("ogs");
#else
        console = spdlog::stdout_color_st("ogs");
#endif
#endif  // defined(_WIN32) && defined(OGS_BUILD_WHEEL)
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
