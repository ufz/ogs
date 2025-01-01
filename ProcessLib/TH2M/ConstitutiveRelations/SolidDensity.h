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

struct SolidDensityModel
{
    void eval(SpaceTimeData const& x_t,
              MediaData const& media_data,
              TemperatureData const& T_data,
              SolidDensityData& solid_density_data) const;

    void dEval(SpaceTimeData const& x_t,
               MediaData const& media_data,
               TemperatureData const& T_data,
               SolidDensityDerivativeData& solid_density_d_data) const;
};

template <int DisplacementDim>
struct SolidDensityModelNonConstantSolidPhaseVolumeFraction
{
    void eval(
        SpaceTimeData const& x_t,
        MediaData const& media_data,
        TemperatureData const& T_data,
        BiotData const& biot,
        StrainData<DisplacementDim> const& strain_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        SolidDensityData& solid_density_data) const;

    void dEval(
        SpaceTimeData const& x_t,
        MediaData const& media_data,
        TemperatureData const& T_data,
        BiotData const& biot,
        StrainData<DisplacementDim> const& strain_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        SolidDensityDerivativeData& solid_density_d_data) const;
};

extern template struct SolidDensityModelNonConstantSolidPhaseVolumeFraction<2>;
extern template struct SolidDensityModelNonConstantSolidPhaseVolumeFraction<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
