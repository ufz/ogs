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
#include "LiquidViscosity.h"
#include "PermeabilityData.h"
#include "ThermoOsmosis.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct DarcyLawData
{
    Eigen::Vector<double, DisplacementDim> v_darcy;

    static auto reflect()
    {
        using Self = DarcyLawData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("velocity",
                                                       &Self::v_darcy);
    }
};

template <int DisplacementDim>
struct DarcyLawModel
{
    explicit DarcyLawModel(
        Eigen::Vector<double, DisplacementDim> const& specific_body_force)
        : b_(specific_body_force)
    {
    }

    void eval(CapillaryPressureData<DisplacementDim> const& p_cap_data,
              LiquidDensityData const& rho_L_data,
              LiquidViscosityData const& mu_L_data,
              PermeabilityData<DisplacementDim> const& perm_data,
              ThermoOsmosisData<DisplacementDim> const& th_osmosis_data,
              DarcyLawData<DisplacementDim>& out) const;

private:
    /// Gravity vector (specific body force).
    Eigen::Vector<double, DisplacementDim> const b_;
};

extern template struct DarcyLawModel<2>;
extern template struct DarcyLawModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
