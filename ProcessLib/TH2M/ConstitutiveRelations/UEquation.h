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
#include "Biot.h"
#include "Bishops.h"
#include "Saturation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct FU2KUpCData
{
    double m = nan;
};

struct FU2KUpCDerivativeData
{
    double dp_cap = nan;
};

struct FU2KUpCModel
{
    void eval(BiotData const& biot_data,
              BishopsData const& chi_S_L,
              FU2KUpCData& fu_2_KupC) const;

    void dEval(BiotData const& biot_data,
               BishopsData const& chi_S_L,
               CapillaryPressureData const& p_cap,
               SaturationDataDeriv const& dS_L_dp_cap,
               FU2KUpCDerivativeData& dfu_2_KupC) const;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
