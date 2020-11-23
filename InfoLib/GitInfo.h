/**
 * \brief  Git information.
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

#include "gitinfolib_export.h"

namespace GitInfoLib
{

namespace GitInfo
{
    extern GITINFOLIB_EXPORT std::string const OGS_VERSION;
    extern GITINFOLIB_EXPORT const std::string git_version_sha1_short;
    extern GITINFOLIB_EXPORT const std::string ogs_version;
}  // namespace
}  // namespace
