// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CouplingScheme.h"

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

namespace ProcessLib
{
bool useMonolithicScheme(std::optional<std::string> const& coupling_scheme,
                         std::string_view const process_name)
{
    if (!coupling_scheme || *coupling_scheme == "monolithic")
    {
        return true;
    }

    if (*coupling_scheme == "monolithic_scheme")
    {
        WARN(
            "The coupling scheme 'monolithic_scheme' of the {:s} process is "
            "deprecated. Use 'monolithic'.",
            process_name);
        return true;
    }

    if (*coupling_scheme == "staggered")
    {
        return false;
    }

    OGS_FATAL(
        "The coupling scheme '{:s}' of the {:s} process is unknown. Use either "
        "'monolithic' or 'staggered'.",
        *coupling_scheme, process_name);
}
}  // namespace ProcessLib
