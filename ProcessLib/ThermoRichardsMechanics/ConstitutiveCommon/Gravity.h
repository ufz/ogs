/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LiquidDensity.h"
#include "Saturation.h"
#include "SolidDensity.h"

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

private:
    Eigen::Vector<double, DisplacementDim> const specific_body_force_;
};

extern template struct GravityModel<2>;
extern template struct GravityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
