// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
/// Returns the number of threads to use for parallel assembly.
/// @return The number of threads set in OGS_ASM_THREADS environment variable or
/// 1 if nothing is set.
int getNumberOfAssemblyThreads();
}  // namespace BaseLib