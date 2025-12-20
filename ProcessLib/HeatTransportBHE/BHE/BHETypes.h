// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <variant>
#include "BHE_1P.h"
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
using BHETypes = std::variant<BHE_1U, BHE_CXA, BHE_CXC, BHE_2U, BHE_1P>;
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
