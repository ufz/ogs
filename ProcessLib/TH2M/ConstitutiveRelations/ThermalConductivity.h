// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "Porosity.h"
#include "Saturation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

template <int DisplacementDim>
struct ThermalConductivityData
{
    GlobalDimMatrix<DisplacementDim> lambda;
};

template <int DisplacementDim>
struct ThermalConductivityDerivativeData
{
    // Currently unused, but there is a comment in TH2MFEM-impl.h referring to
    // this matrix
    // GlobalDimMatrix<DisplacementDim> dlambda_dp_GR;
    GlobalDimMatrix<DisplacementDim> dlambda_dp_cap;
    GlobalDimMatrix<DisplacementDim> dlambda_dT;
};

template <int DisplacementDim>
struct ThermalConductivityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              TemperatureData const& T_data, PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              ThermalConductivityData<DisplacementDim>&
                  thermal_conductivity_data) const;

    void dEval(SpaceTimeData const& x_t, MediaData const& media_data,
               TemperatureData const& T_data, PorosityData const& porosity_data,
               PorosityDerivativeData const& porosity_d_data,
               SaturationData const& S_L_data,
               ThermalConductivityDerivativeData<DisplacementDim>&
                   thermal_conductivity_d_data) const;
};

extern template struct ThermalConductivityModel<2>;
extern template struct ThermalConductivityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
