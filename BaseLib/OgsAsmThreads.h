// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
/// Reads the integer value from the given environment variable and returns it.
/// Returns 1 if the variable is unset. Calls OGS_FATAL if the variable is set
/// but empty, cannot be parsed as an integer, or its value is less than 1.
int getNumberOfThreadsFromEnv(char const* env_var_name);

/// Returns the number of threads to use for parallel assembly.
/// @return The number of threads set in OGS_ASM_THREADS environment variable or
/// 1 if nothing is set.
int getNumberOfAssemblyThreads();

/// Returns the maximum of getNumberOfAssemblyThreads() and
/// omp_get_max_threads().
int getNumberOfThreads();
}  // namespace BaseLib
