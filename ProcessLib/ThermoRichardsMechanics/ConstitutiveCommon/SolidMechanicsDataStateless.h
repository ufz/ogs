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

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SolidMechanicsDataStateless
{
    KelvinMatrix<DisplacementDim> stiffness_tensor =
        KV::KMnan<DisplacementDim>();
    KelvinVector<DisplacementDim> J_uT_BT_K_N = KV::KVnan<DisplacementDim>();
    KelvinVector<DisplacementDim> J_up_BT_K_N = KV::KVnan<DisplacementDim>();
};
}  // namespace ProcessLib::ThermoRichardsMechanics
