// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "BaseLib/Logging.h"
#include "LisWrapper.h"

namespace MathLib
{

/**
 * check Lis error codes
 *
 * \param err   Lis error code
 * \return success or not
 */
inline bool checkLisError(int err)
{
    bool ok = (err == LIS_SUCCESS);
    if (!ok)
    {
        ERR("***ERROR: Lis error code = {:d}", err);
    }
    return ok;
}

}  // namespace MathLib
