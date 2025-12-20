// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "Biot.h"
#include "Bishops.h"
#include "Saturation.h"
#include "SolidMechanics.h"
#include "SolidThermalExpansion.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct FU1KUTDerivativeData
{
    KelvinVector<DisplacementDim> dT;
};

template <int DisplacementDim>
struct FU1KUTModel
{
    void dEval(
        SolidMechanicsDataStateless<DisplacementDim> const& s_mech_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        FU1KUTDerivativeData<DisplacementDim>& dfu_1_KuT) const;
};

extern template struct FU1KUTModel<2>;
extern template struct FU1KUTModel<3>;

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
