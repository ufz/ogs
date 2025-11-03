/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace BaseLib
{
/// Returns the number of threads to use for parallel assembly.
/// @return The number of threads set in OGS_ASM_THREADS environment variable or
/// 1 if nothing is set.
int getNumberOfAssemblyThreads();
}  // namespace BaseLib