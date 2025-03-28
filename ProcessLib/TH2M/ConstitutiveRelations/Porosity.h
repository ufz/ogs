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
#include "Biot.h"
#include "Bishops.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/Reflection/ReflectionData.h"
#include "Saturation.h"
#include "SolidCompressibility.h"
#include "SolidThermalExpansion.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct PorosityDerivativeData
{
    double dphi_dT = nan;
    // dphi_G_dp_GR = -ds_L_dp_GR * phi = 0;
    // dphi_L_dp_GR = ds_L_dp_GR * phi = 0;
    double dphi_L_dp_cap = nan;
    // dphi_G_dp_cap = -dphi_L_dp_cap
};

struct PorosityData
{
    double phi = nan;

    static auto reflect()
    {
        using Self = PorosityData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{R::makeReflectionData("porosity", &Self::phi)};
    }
};

template <int DisplacementDim>
struct PorosityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
        SaturationData const& S_L_data,
        PrevState<SaturationData> const& S_L_prev_data,
        CapillaryPressureData const& p_cap, GasPressureData const& p_GR,
        BishopsData const& chi_S_L,
        PrevState<BishopsData> const& chi_S_L_prev, BiotData const& biot,
        SolidCompressibilityData const& solid_compressibility,
        StrainData<DisplacementDim> const& eps_data,
        PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
        PrevState<PorosityData> const& porosity_prev_data,
        PorosityData& porosity_data) const;

    void dEval(
        SpaceTimeData const& x_t, MediaData const& media_data,
        PorosityData const& porosity_data,
        SaturationDataDeriv const& dS_L_dp_cap,
        PorosityDerivativeData& porosity_d_data) const;
};



extern template struct PorosityModel<2>;
extern template struct PorosityModel<3>;

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
