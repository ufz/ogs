/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Base.h"

namespace ProcessLib::RichardsMechanics
{
struct SaturationSecantDerivative
{
    double DeltaS_L_Deltap_cap = nan;
};
}  // namespace ProcessLib::RichardsMechanics
