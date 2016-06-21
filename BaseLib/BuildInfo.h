/**
 * \brief  Build information.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BASE_LIB_BUILD_INFO_H_
#define BASE_LIB_BUILD_INFO_H_

#include <string>
#ifdef _MSC_VER
#include "baselib_export.h"
#else
#define BASELIB_EXPORT
#endif

namespace BaseLib
{

namespace BuildInfo
{
    extern const std::string build_timestamp;

    extern const std::string cmake_cxx_compiler;
    extern const std::string cmake_cxx_flags;
    extern const std::string cmake_cxx_flags_release;
    extern const std::string cmake_cxx_flags_debug;

    extern const std::string git_version_sha1;
    extern const std::string git_version_sha1_short;

	extern const BASELIB_EXPORT std::string git_describe;
    extern const std::string ogs_version;

    extern const std::string source_path;
    extern const BASELIB_EXPORT std::string data_path;
    extern const std::string data_binary_path;
    extern const BASELIB_EXPORT std::string tests_tmp_path;
}
}

#endif  // BASE_LIB_BUILD_INFO_H_
