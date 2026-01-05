// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "LiquidDensity.h"
#include "Saturation.h"
#include "SolidDensity.h"
#include "SpecificBodyForceData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct GravityData
{
    GlobalDimVector<DisplacementDim> volumetric_body_force;
    GlobalDimVector<DisplacementDim> J_up_HT_V_N;
};

template <int DisplacementDim>
struct GravityModel
{
    explicit GravityModel(
        Eigen::Vector<double, DisplacementDim> const& specific_body_force)
        : specific_body_force_(specific_body_force)
    {
    }

    void eval(PorosityData const& poro_data,
              SolidDensityData const& rho_S_data,
              LiquidDensityData const& rho_L_data,
              SaturationData const& S_L_data,
              SaturationDataDeriv const& dS_L_data,
              GravityData<DisplacementDim>& out) const;

    static GravityModel create(
        SpecificBodyForceData<DisplacementDim> const& specific_body_force_data)
    {
        return GravityModel{specific_body_force_data.specific_body_force};
    }

private:
    Eigen::Vector<double, DisplacementDim> const specific_body_force_;
};

extern template struct GravityModel<2>;
extern template struct GravityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
