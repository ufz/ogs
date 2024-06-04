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

#include "Density.h"
#include "LiquidDensity.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Biot.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/LiquidViscosity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/PermeabilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Porosity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Saturation.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidCompressibilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStress_StrainTemperature/SolidMechanics.h"
#include "SaturationSecantDerivative.h"
#include "StiffnessTensor.h"

namespace ProcessLib::RichardsMechanics
{
/// Data whose state must be tracked by the TRM process.
template <int DisplacementDim>
using StatefulData = std::tuple<
    StrainData<DisplacementDim>,
    ProcessLib::ThermoRichardsMechanics::ConstitutiveStress_StrainTemperature::
        EffectiveStressData<DisplacementDim>>;

template <int DisplacementDim>
using StatefulDataPrev = ProcessLib::ConstitutiveRelations::PrevStateOf<
    StatefulData<DisplacementDim>>;

/// Data that is needed for output purposes, but not directly for the assembly.
template <int DisplacementDim>
using OutputData = std::tuple<>;

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
using ConstitutiveData = std::tuple<
    // TODO (CL) check if all that data should stay here
    StiffnessTensor<DisplacementDim>,
    ProcessLib::ThermoRichardsMechanics::PorosityData, Density, LiquidDensity,
    ProcessLib::ThermoRichardsMechanics::BiotData,
    ProcessLib::ThermoRichardsMechanics::SaturationDataDeriv,
    ProcessLib::ThermoRichardsMechanics::LiquidViscosityData,
    ProcessLib::ThermoRichardsMechanics::SolidCompressibilityData,
    ProcessLib::ThermoRichardsMechanics::BishopsData,
    PrevState<ProcessLib::ThermoRichardsMechanics::BishopsData>,
    ProcessLib::ThermoRichardsMechanics::PermeabilityData<DisplacementDim>,
    SaturationSecantDerivative>;

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
using ConstitutiveTempData = std::tuple<>;
}  // namespace ProcessLib::RichardsMechanics
