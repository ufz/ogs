// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "FluidDensity.h"
#include "PermeabilityData.h"
#include "ProcessLib/Reflection/ReflectionData.h"
#include "Viscosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct DarcyVelocityData
{
    GlobalDimVector<DisplacementDim> w_GS;
    GlobalDimVector<DisplacementDim> w_LS;

    static auto reflect()
    {
        using Self = DarcyVelocityData<DisplacementDim>;
        namespace R = ProcessLib::Reflection;

        return std::tuple{
            R::makeReflectionData("velocity_gas", &Self::w_GS),
            R::makeReflectionData("velocity_liquid", &Self::w_LS)};
    }
};

template <int DisplacementDim>
struct DarcyVelocityModel
{
    void eval(CapillaryPressureGradientData<DisplacementDim> const& grad_p_cap,
              FluidDensityData const& fluid_density_data,
              GasPressureGradientData<DisplacementDim> const& grad_p_GR,
              PermeabilityData<DisplacementDim> const& permeability_data,
              SpecificBodyForceData<DisplacementDim> const& specific_body_force,
              ViscosityData const& viscosity_data,
              DarcyVelocityData<DisplacementDim>& darcy_velocity_data) const;
};

extern template struct DarcyVelocityModel<2>;
extern template struct DarcyVelocityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
