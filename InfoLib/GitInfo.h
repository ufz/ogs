// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
}  // namespace GitInfo
}  // namespace GitInfoLib
