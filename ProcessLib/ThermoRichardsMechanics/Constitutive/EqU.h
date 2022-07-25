/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LiquidDensity.h"
#include "Porosity.h"
#include "Saturation.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct EqUData
{
    GlobalDimVector<DisplacementDim> J_up_HT_V_N = DVnan<DisplacementDim>();
    double J_up_X_BTI2N = nan;
};

template <int DisplacementDim>
struct EqUModel
{
    explicit EqUModel(
        Eigen::Vector<double, DisplacementDim> const& specific_body_force)
        : b_(specific_body_force)
    {
    }

    void eval(CapillaryPressureData<DisplacementDim> const& p_cap_data,
              SaturationDataDeriv const& dS_L_data,
              BiotData const& biot_data,
              BishopsData const& bishops_data,
              LiquidDensityData const& rho_L_data,
              PorosityData const& poro_data,
              EqUData<DisplacementDim>& out) const;

private:
    /// Gravity vector (specific body force).
    Eigen::Vector<double, DisplacementDim> const b_;
};

extern template struct EqUModel<2>;
extern template struct EqUModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
