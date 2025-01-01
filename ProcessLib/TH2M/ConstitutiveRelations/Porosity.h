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
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/Reflection/ReflectionData.h"
#include "Saturation.h"
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

struct PorosityModel
{
    void eval(SpaceTimeData const& x_t,
              MediaData const& media_data,
              PorosityData& porosity_data) const;

    void dEval(SpaceTimeData const& x_t,
               MediaData const& media_data,
               PorosityData const& porosity_data,
               SaturationDataDeriv const& dS_L_dp_cap,
               PorosityDerivativeData& porosity_d_data) const;
};

template <int DisplacementDim>
struct PorosityModelNonConstantSolidPhaseVolumeFraction
{
    void eval(
        SpaceTimeData const& x_t,
        MediaData const& media_data,
        BiotData const& biot,
        StrainData<DisplacementDim> const& strain_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        PorosityData& porosity_data) const;

    void dEval(
        SpaceTimeData const& x_t,
        MediaData const& media_data,
        PorosityData const& porosity_data,
        SaturationDataDeriv const& dS_L_dp_cap,
        BiotData const& biot,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        StrainData<DisplacementDim> const& strain_data,
        PorosityDerivativeData& porosity_d_data) const;
};

extern template struct PorosityModelNonConstantSolidPhaseVolumeFraction<2>;
extern template struct PorosityModelNonConstantSolidPhaseVolumeFraction<3>;

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
