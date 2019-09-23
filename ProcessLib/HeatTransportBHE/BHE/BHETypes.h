/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <variant>
#include "BHE_1U.h"
#include "BHE_2U.h"
#include "BHE_CXA.h"
#include "BHE_CXC.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
using BHETypes = std::variant<BHE_1U, BHE_CXA, BHE_CXC, BHE_2U>;
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
