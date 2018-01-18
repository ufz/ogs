/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessOutput.h"

namespace ProcessLib
{
//! Creates an instance of ProcessOutput from config.
ProcessOutput createProcessOutput(BaseLib::ConfigTree const& output_config);

}  // namespace ProcessLib
