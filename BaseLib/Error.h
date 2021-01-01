/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstdlib>
#include "Logging.h"

#ifdef OGS_FATAL_ABORT
#define OGS_FATAL(...)                                                      \
    {                                                                       \
        BaseLib::console->critical("{}:{} {}() ", __FILE__, __LINE__,       \
                                   __FUNCTION__, fmt::format(__VA_ARGS__)); \
        std::abort();                                                       \
    }
#else  // OGS_FATAL_ABORT
#include <stdexcept>
#define OGS_FATAL(...)                                                      \
    {                                                                       \
        BaseLib::console->critical("{}:{} {}() ", __FILE__, __LINE__,       \
                                   __FUNCTION__, fmt::format(__VA_ARGS__)); \
        throw std::runtime_error(fmt::format(__VA_ARGS__));                 \
    }
#endif  // OGS_FATAL_ABORT

