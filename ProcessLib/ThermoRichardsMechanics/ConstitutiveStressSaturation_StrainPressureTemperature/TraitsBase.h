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
