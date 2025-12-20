// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <tclap/CmdLine.h>

#include <string>

namespace BaseLib
{

TCLAP::ValueArg<std::string> makeLogLevelArg();

}  // namespace BaseLib
