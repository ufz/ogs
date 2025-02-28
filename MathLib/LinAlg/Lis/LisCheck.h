/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <lis.h>

#include <vector>

#include "BaseLib/Logging.h"

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
