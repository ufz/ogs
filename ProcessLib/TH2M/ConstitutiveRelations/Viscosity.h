/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
