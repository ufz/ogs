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

#include "ProcessLib/Reflection/ReflectionData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/DarcyLaw.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqP.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqT.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EquivalentPlasticStrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Gravity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Porosity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Saturation.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidDensity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Swelling.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMHeatStorageAndFlux.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMStorage.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMVaporDiffusion.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TotalStressData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TransportPorosity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStressSaturation_StrainPressureTemperature/SolidMechanics.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
/// Data whose state must be tracked by the TRM process.
template <int DisplacementDim>
using StatefulData =
    std::tuple<SaturationData, PorosityData, TransportPorosityData,
               StrainData<DisplacementDim>,
               // TODO swelling is not used in this constitutive setting, but
               // the stateful data is set in setInitialConditionsConcrete()
               // Once this has been refactored, the swelling code can be moved
               // from the common to the original subdirectory for constitutive
               // settings
               SwellingDataStateful<DisplacementDim>,
               SolidMechanicsDataStateful<DisplacementDim>,
               TotalStressData<DisplacementDim>>;

// TODO convert to tuple
/// Data whose state must be tracked by the TRM process.
template <int DisplacementDim>
struct StatefulDataPrev
{
    PrevState<SaturationData> S_L_data;
    PrevState<PorosityData> poro_data;
    PrevState<TransportPorosityData> transport_poro_data;
    PrevState<StrainData<DisplacementDim>> eps_data;
    PrevState<SwellingDataStateful<DisplacementDim>> swelling_data;
    PrevState<SolidMechanicsDataStateful<DisplacementDim>> s_mech_data;
    PrevState<TotalStressData<DisplacementDim>> total_stress_data;

    StatefulDataPrev<DisplacementDim>& operator=(
        StatefulData<DisplacementDim> const& state)
    {
        S_L_data = std::get<SaturationData>(state);
        poro_data = std::get<PorosityData>(state);
        transport_poro_data = std::get<TransportPorosityData>(state);
        eps_data = std::get<StrainData<DisplacementDim>>(state);
        swelling_data = std::get<SwellingDataStateful<DisplacementDim>>(state);
        s_mech_data =
            std::get<SolidMechanicsDataStateful<DisplacementDim>>(state);
        total_stress_data = std::get<TotalStressData<DisplacementDim>>(state);

        return *this;
    }
};

/// Data that is needed for output purposes, but not directly for the assembly.
template <int DisplacementDim>
using OutputData = std::tuple<DarcyLawData<DisplacementDim>,
                              LiquidDensityData,
                              LiquidViscosityData,
                              SolidDensityData,
                              PermeabilityData<DisplacementDim>>;

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
using ConstitutiveData =
    std::tuple<SolidMechanicsDataStateless<DisplacementDim>,
               GravityData<DisplacementDim>,
               TRMHeatStorageAndFluxData<DisplacementDim>,
               TRMVaporDiffusionData<DisplacementDim>,
               TRMStorageData,
               EqPData<DisplacementDim>,
               EqTData<DisplacementDim>,
               ThermoOsmosisData<DisplacementDim>>;

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
using ConstitutiveTempData =
    std::tuple<ElasticTangentStiffnessData<DisplacementDim>,
               BiotData,
               SolidCompressibilityData,
               SaturationDataDeriv,
               BishopsData,
               // TODO why not usual state tracking for that?
               PrevState<BishopsData>,
               SolidThermalExpansionData<DisplacementDim>,
               FluidThermalExpansionData,
               EquivalentPlasticStrainData>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
