/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using MicroSaturation = BaseLib::StrongType<double, struct MicroSaturationTag>;
}  // namespace ProcessLib::RichardsMechanics
