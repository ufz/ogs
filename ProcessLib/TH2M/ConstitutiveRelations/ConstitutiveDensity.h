// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
/// Real constituent partial densities.
struct ConstituentDensityData
{
    // Gas phase density
    double rho_C_GR = nan;
    double rho_W_GR = nan;

    // Liquid phase density
    double rho_C_LR = nan;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
