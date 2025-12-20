// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "OgsAsmThreads.h"

#include <cstdlib>
#include <sstream>

#include "Error.h"
#include "StringTools.h"

namespace BaseLib
{
int getNumberOfAssemblyThreads()
{
    char const* const num_threads_env = std::getenv("OGS_ASM_THREADS");

    if (!num_threads_env)
    {
        return 1;
    }

    if (std::strlen(num_threads_env) == 0)
    {
        OGS_FATAL("The environment variable OGS_ASM_THREADS is set but empty.");
    }

    std::string num_threads_str{num_threads_env};
    BaseLib::trim(num_threads_str);

    std::istringstream num_threads_iss{num_threads_str};
    int num_threads = -1;

    num_threads_iss >> num_threads;

    if (!num_threads_iss)
    {
        OGS_FATAL("Error parsing OGS_ASM_THREADS (= \"{}\").", num_threads_env);
    }

    if (!num_threads_iss.eof())
    {
        OGS_FATAL(
            "Error parsing OGS_ASM_THREADS (= \"{}\"): not read entirely, the "
            "remainder is \"{}\"",
            num_threads_env,
            num_threads_iss.str().substr(num_threads_iss.tellg()));
    }

    if (num_threads < 1)
    {
        OGS_FATAL(
            "You asked (via OGS_ASM_THREADS) to assemble with {} < 1 thread.",
            num_threads);
    }

    return num_threads;
}
}  // namespace BaseLib