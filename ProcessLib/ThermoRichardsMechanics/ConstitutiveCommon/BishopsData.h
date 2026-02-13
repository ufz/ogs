// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"  // for nan

namespace ProcessLib::ThermoRichardsMechanics
{
struct BishopsData
{
    double chi_S_L = nan;
    double dchi_dS_L = nan;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
