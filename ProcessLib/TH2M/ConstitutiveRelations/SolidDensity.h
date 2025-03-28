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
#include "Porosity.h"
#include "ProcessLib/ConstitutiveRelations/EffectiveStressData.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/Reflection/ReflectionData.h"
#include "SolidThermalExpansion.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct SolidDensityDerivativeData
{
    double drho_SR_dT = nan;
};

struct SolidDensityData
{
    double rho_SR = nan;

    static auto reflect()
    {
        using Self = SolidDensityData;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("solid_density", &Self::rho_SR)};
    }
};

template <int DisplacementDim>
struct SolidDensityModel
{
    void eval(SpaceTimeData const& x_t,
        MediaData const& media_data,
        TemperatureData const& T_data,
        ProcessLib::ConstitutiveRelations::EffectiveStressData<
               DisplacementDim> const& sigma_eff_data,
        CapillaryPressureData const& p_cap,
        GasPressureData const& p_GR,
        BishopsData const& chi_S_L,
        PorosityData const& poro_data,
        SolidDensityData& solid_density_data) const;

    void dEval(SpaceTimeData const& x_t,
        MediaData const& media_data,
        TemperatureData const& T_data,
        ProcessLib::ConstitutiveRelations::EffectiveStressData<
                                        DisplacementDim> const& sigma_eff_data,
        CapillaryPressureData const& p_cap,
        GasPressureData const& p_GR,
        BishopsData const& chi_S_L,
        PorosityData const& poro_data,
        SolidDensityDerivativeData& solid_density_d_data) const;
};

extern template struct SolidDensityModel<2>;
extern template struct SolidDensityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
