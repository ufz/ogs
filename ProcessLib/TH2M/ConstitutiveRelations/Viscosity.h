// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "MassMoleFractions.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct ViscosityData
{
    double mu_GR = nan;
    double mu_LR = nan;
};

struct ViscosityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              TemperatureData const& T_data,
              MassMoleFractionsData const& mass_mole_fractions_data,
              ViscosityData& viscosity_data) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
