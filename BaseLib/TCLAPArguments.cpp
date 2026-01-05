// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TCLAPArguments.h"

#include <vector>

namespace BaseLib
{

TCLAP::ValueArg<std::string> makeLogLevelArg()
{
    static std::vector<std::string> allowed_log_levels{"none", "error", "warn",
                                                       "info", "debug", "all"};
    static auto* allowed_log_levels_vals =
        new TCLAP::ValuesConstraint<std::string>(allowed_log_levels);

    return TCLAP::ValueArg<std::string>(
        "l", "log-level", "the verbosity of logging messages", false,
#ifdef NDEBUG
        "info",
#else
        "all",
#endif
        allowed_log_levels_vals);
}

}  // namespace BaseLib
