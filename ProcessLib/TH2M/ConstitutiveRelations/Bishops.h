/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "Saturation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct BishopsData
{
    double chi_S_L = nan;
    double dchi_dS_L = nan;
};

struct BishopsModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SaturationData const& S_L_data, BishopsData& out) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
