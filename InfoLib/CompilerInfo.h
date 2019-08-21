/**
 * \brief  Build information.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

#include "compilerinfolib_export.h"

namespace CompilerInfoLib
{

namespace CompilerInfo
{
    extern COMPILERINFOLIB_EXPORT const std::string cmake_cxx_compiler; // all not used
    extern COMPILERINFOLIB_EXPORT const std::string cmake_cxx_flags;
    extern COMPILERINFOLIB_EXPORT const std::string cmake_cxx_flags_release;
    extern COMPILERINFOLIB_EXPORT const std::string cmake_cxx_flags_debug;
    }  // namespace BuildInfo
    }  // namespace BaseLib
