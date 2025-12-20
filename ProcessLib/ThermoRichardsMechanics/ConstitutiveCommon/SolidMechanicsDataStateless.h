// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
