/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <tclap/CmdLine.h>

#include <string>

namespace BaseLib
{

TCLAP::ValueArg<std::string> makeLogLevelArg();

}  // namespace BaseLib
