/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifdef OGS_FATAL_ABORT

#include <cstdlib>
#include <logog/include/logog.hpp>

namespace BaseLib
{
namespace detail
{
template <typename Msg>
[[noreturn]] bool error_impl(Msg&& msg)
{
    ERR("%s", msg.data());
    std::abort();
}

}  // namespace detail

}  // namespace BaseLib

#else  // OGS_FATAL_ABORT

#include <stdexcept>

namespace BaseLib
{
namespace detail
{
template <typename Msg>
[[noreturn]] bool error_impl(Msg&& msg)
{
    throw std::runtime_error(std::forward<Msg>(msg));
}

}  // namespace detail

}  // namespace BaseLib

#endif  // OGS_FATAL_ABORT

#include "FileTools.h"
#include "StringTools.h"

#define OGS_STR(x) #x
#define OGS_STRINGIFY(x) OGS_STR(x)
#define OGS_LOCATION                              \
    " at " + BaseLib::extractBaseName(__FILE__) + \
        ", line " OGS_STRINGIFY(__LINE__)

#define OGS_FATAL(fmt, ...)                                           \
    BaseLib::detail::error_impl(BaseLib::format(fmt, ##__VA_ARGS__) + \
                                OGS_LOCATION)
