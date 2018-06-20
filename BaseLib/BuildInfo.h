/**
 * \brief  Build information.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

#include "baselib_export.h"

namespace BaseLib
{

namespace BuildInfo
{
    extern BASELIB_EXPORT const std::string build_timestamp;

    extern BASELIB_EXPORT const std::string cmake_cxx_compiler;
    extern BASELIB_EXPORT const std::string cmake_cxx_flags;
    extern BASELIB_EXPORT const std::string cmake_cxx_flags_release;
    extern BASELIB_EXPORT const std::string cmake_cxx_flags_debug;

    extern BASELIB_EXPORT const std::string git_version_sha1;
    extern BASELIB_EXPORT const std::string git_version_sha1_short;

    extern BASELIB_EXPORT const std::string git_describe;
    extern BASELIB_EXPORT const std::string ogs_version;

    extern BASELIB_EXPORT const std::string source_path;
    extern BASELIB_EXPORT const std::string data_path;
    extern BASELIB_EXPORT const std::string data_binary_path;
    extern BASELIB_EXPORT const std::string tests_tmp_path;
}
}
