/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/variant.hpp>
#include "BHE_1U.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
using BHETypes = boost::variant<BHE_1U>;
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
