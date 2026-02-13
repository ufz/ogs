// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CapillaryPressureData.h"
#include "EqPData.h"
#include "EqTData.h"
#include "GravityData.h"
#include "PermeabilityData.h"
#include "SolidThermalExpansionData.h"
#include "TRMHeatStorageAndFluxData.h"
#include "TRMVaporDiffusionData.h"
#include "TemperatureData.h"
#include "ThermoOsmosisData.h"
#include "TotalStressData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
// Explicit instantiations for templated Data structs used in 2D and 3D builds.
template struct SolidThermalExpansionData<2>;
template struct SolidThermalExpansionData<3>;

template struct CapillaryPressureData<2>;
template struct CapillaryPressureData<3>;

template struct TotalStressData<2>;
template struct TotalStressData<3>;

template struct TRMVaporDiffusionData<2>;
template struct TRMVaporDiffusionData<3>;

template struct TRMHeatStorageAndFluxData<2>;
template struct TRMHeatStorageAndFluxData<3>;

template struct ThermoOsmosisData<2>;
template struct ThermoOsmosisData<3>;

template struct EqPData<2>;
template struct EqPData<3>;

template struct EqTData<2>;
template struct EqTData<3>;

template struct PermeabilityData<2>;
template struct PermeabilityData<3>;

template struct TemperatureData<2>;
template struct TemperatureData<3>;

template struct GravityData<2>;
template struct GravityData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
