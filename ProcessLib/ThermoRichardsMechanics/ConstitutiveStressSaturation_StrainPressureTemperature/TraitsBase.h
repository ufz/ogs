// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/SolidModels/MFront/MFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
namespace MSM = MaterialLib::Solids::MFront;

template <int DisplacementDim>
using SolidConstitutiveRelation =
    MSM::MFrontGeneric<DisplacementDim,
                       boost::mp11::mp_list<MSM::Strain, MSM::LiquidPressure>,
                       boost::mp11::mp_list<MSM::Stress, MSM::Saturation>,
                       boost::mp11::mp_list<MSM::Temperature>>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
