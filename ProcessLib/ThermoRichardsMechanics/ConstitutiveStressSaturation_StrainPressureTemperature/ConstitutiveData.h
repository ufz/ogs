// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <tuple>

#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/ConstitutiveRelations/EffectiveStressData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Biot.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/BishopsData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/DarcyLawData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqPData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqTData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EquivalentPlasticStrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/FluidThermalExpansionData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/GravityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/LiquidDensityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/LiquidViscosityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/PermeabilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/PorosityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SaturationData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidCompressibilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidDensityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidMechanicsDataStateless.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidThermalExpansionData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMHeatStorageAndFluxData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMStorageData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMVaporDiffusionData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/ThermoOsmosisData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TotalStressData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TransportPorosityData.h"
#include "SolidMechanics.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
// TODO directly declare these type aliases in Traits.h
/// Data whose state must be tracked by the TRM process.
template <int DisplacementDim>
using StatefulData =
    std::tuple<SaturationData, PorosityData, TransportPorosityData,
               StrainData<DisplacementDim>, TotalStressData<DisplacementDim>>;

template <int DisplacementDim>
using StatefulDataPrev = PrevStateOf<StatefulData<DisplacementDim>>;

/// Data that is needed for output purposes, but not directly for the assembly.
template <int DisplacementDim>
using OutputData = std::tuple<DarcyLawData<DisplacementDim>, LiquidDensityData,
                              LiquidViscosityData, SolidDensityData,
                              PermeabilityData<DisplacementDim>>;

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
using ConstitutiveData =
    std::tuple<SolidMechanicsDataStateless<DisplacementDim>,
               EffectiveStressData<DisplacementDim>,
               GravityData<DisplacementDim>,
               TRMHeatStorageAndFluxData<DisplacementDim>,
               TRMVaporDiffusionData<DisplacementDim>, TRMStorageData,
               EqPData<DisplacementDim>, EqTData<DisplacementDim>,
               ThermoOsmosisData<DisplacementDim>>;

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
using ConstitutiveTempData = std::tuple<
    BiotData, SolidCompressibilityData, SaturationDataDeriv, BishopsData,
    // TODO why not usual state tracking for that?
    PrevState<BishopsData>, SolidThermalExpansionData<DisplacementDim>,
    FluidThermalExpansionData, EquivalentPlasticStrainData>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
