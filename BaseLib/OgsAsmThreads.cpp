// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "OgsAsmThreads.h"

#include <algorithm>
#include <cstdlib>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Error.h"
#include "StringTools.h"

namespace BaseLib
{
int getNumberOfThreadsFromEnv(char const* env_var_name)
{
    char const* const num_threads_env = std::getenv(env_var_name);

    if (!num_threads_env)
    {
        return 1;
    }

    if (std::strlen(num_threads_env) == 0)
    {
        OGS_FATAL("The environment variable {} is set but empty.",
                  env_var_name);
    }

    std::string num_threads_str{num_threads_env};
    BaseLib::trim(num_threads_str);

    std::istringstream num_threads_iss{num_threads_str};
    int num_threads = -1;

    num_threads_iss >> num_threads;

    if (!num_threads_iss)
    {
        OGS_FATAL("Error parsing {} (= \"{}\").", env_var_name,
                  num_threads_env);
    }

    if (!num_threads_iss.eof())
    {
        OGS_FATAL(
            "Error parsing {} (= \"{}\"): not read entirely, the "
            "remainder is \"{}\"",
            env_var_name,
            num_threads_env,
            num_threads_iss.str().substr(num_threads_iss.tellg()));
    }

    if (num_threads < 1)
    {
        OGS_FATAL("You asked (via {}) to use {} < 1 thread.", env_var_name,
                  num_threads);
    }

    return num_threads;
}

int getNumberOfAssemblyThreads()
{
    return getNumberOfThreadsFromEnv("OGS_ASM_THREADS");
}

int getNumberOfThreads()
{
#ifdef _OPENMP
    int const num_omp_threads = omp_get_max_threads();
#else
    int const num_omp_threads = 1;
#endif
    return std::max(getNumberOfAssemblyThreads(), num_omp_threads);
}
}  // namespace BaseLib
