// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "OgsChemThreads.h"

#include <cstdlib>

#include "OgsAsmThreads.h"

namespace BaseLib
{
int getNumberOfChemistryThreads()
{
    if (std::getenv("OGS_CHEM_THREADS"))
    {
        return getNumberOfThreadsFromEnv("OGS_CHEM_THREADS");
    }
    return getNumberOfAssemblyThreads();
}
}  // namespace BaseLib
