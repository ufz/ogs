// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdlib>

#include "Logging.h"

#ifdef OGS_FATAL_ABORT
#define OGS_FATAL(...)                                                  \
    {                                                                   \
        auto const formatted_va_args = fmt::format(__VA_ARGS__);        \
        BaseLib::console->critical("{}:{} {}() {}", __FILE__, __LINE__, \
                                   __FUNCTION__, formatted_va_args);    \
        BaseLib::console->flush();                                      \
        std::abort();                                                   \
    }
#else  // OGS_FATAL_ABORT
#include <stdexcept>
#define OGS_FATAL(...)                                                  \
    {                                                                   \
        auto const formatted_va_args = fmt::format(__VA_ARGS__);        \
        BaseLib::console->critical("{}:{} {}() {}", __FILE__, __LINE__, \
                                   __FUNCTION__, formatted_va_args);    \
        BaseLib::console->flush();                                      \
        throw std::runtime_error(formatted_va_args);                    \
    }
#endif  // OGS_FATAL_ABORT
