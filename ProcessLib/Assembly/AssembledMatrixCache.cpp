// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AssembledMatrixCache.h"

namespace ProcessLib
{
AssembledMatrixCache::AssembledMatrixCache(bool const is_linear,
                                           bool const use_monolithic_scheme)
    : is_linear_{is_linear && use_monolithic_scheme}
{
    if (is_linear && !use_monolithic_scheme)
    {
        OGS_FATAL(
            "You requested to assemble only once in combination with staggered "
            "coupling. This use case is not yet implemented.");
    }

    if (is_linear_)
    {
        WARN(
            "You specified that the process simulated by OGS is linear. With "
            "that optimization the process will be assembled only once and the "
            "non-linear solver will do only one iteration per time step. No "
            "non-linearities will be resolved and OGS will not detect if there "
            "are any non-linearities. It is your responsibility to ensure that "
            "the assembled equation systems are linear, indeed! There is no "
            "safety net!");
    }
}
}  // namespace ProcessLib
