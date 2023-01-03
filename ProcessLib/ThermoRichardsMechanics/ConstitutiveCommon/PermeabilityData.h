/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct PermeabilityData
{
    double k_rel;
    double dk_rel_dS_L;
    GlobalDimMatrix<DisplacementDim> Ki_over_mu;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
