/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStress_StrainTemperature/SolidMechanics.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
/// Data whose state must be tracked by the TRM process.
template <int DisplacementDim>
struct StatefulData
{
    SaturationData S_L_data;
    PorosityData poro_data;
    TransportPorosityData transport_poro_data;
    StrainData<DisplacementDim> eps_data;
    SwellingDataStateful<DisplacementDim> swelling_data;
    SolidMechanicsDataStateful<DisplacementDim> s_mech_data;

    static auto reflect()
    {
        using Self = StatefulData<DisplacementDim>;

        return Reflection::reflectWithoutName(&Self::S_L_data,
                                              &Self::poro_data,
                                              &Self::transport_poro_data,
                                              &Self::eps_data,
                                              &Self::swelling_data,
                                              &Self::s_mech_data);
    }
};

template <int DisplacementDim>
struct StatefulDataPrev
{
    PrevState<SaturationData> S_L_data;
    PrevState<PorosityData> poro_data;
    PrevState<TransportPorosityData> transport_poro_data;
    PrevState<StrainData<DisplacementDim>> eps_data;
    PrevState<SwellingDataStateful<DisplacementDim>> swelling_data;
    PrevState<SolidMechanicsDataStateful<DisplacementDim>> s_mech_data;

    StatefulDataPrev<DisplacementDim>& operator=(
        StatefulData<DisplacementDim> const& state)
    {
        S_L_data = state.S_L_data;
        poro_data = state.poro_data;
        transport_poro_data = state.transport_poro_data;
        eps_data = state.eps_data;
        swelling_data = state.swelling_data;
        s_mech_data = state.s_mech_data;

        return *this;
    }
};

/// Data that is needed for output purposes, but not directly for the assembly.
template <int DisplacementDim>
struct OutputData
{
    DarcyLawData<DisplacementDim> darcy_data;
    LiquidDensityData rho_L_data;
    LiquidViscosityData mu_L_data;
    SolidDensityData rho_S_data;

    static auto reflect()
    {
        using Self = OutputData<DisplacementDim>;

        return Reflection::reflectWithoutName(&Self::darcy_data,
                                              &Self::rho_L_data,
                                              &Self::mu_L_data,
                                              &Self::rho_S_data);
    }
};

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
struct ConstitutiveData
{
    SolidMechanicsDataStateless<DisplacementDim> s_mech_data;
    TotalStressData<DisplacementDim> total_stress_data;
    GravityData<DisplacementDim> grav_data;
    TRMHeatStorageAndFluxData<DisplacementDim> heat_data;
    TRMVaporDiffusionData<DisplacementDim> vap_data;
    TRMStorageData storage_data;
    EqPData<DisplacementDim> eq_p_data;
    EqTData<DisplacementDim> eq_T_data;
    ThermoOsmosisData<DisplacementDim> th_osmosis_data;
};

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
struct ConstitutiveTempData
{
    SwellingDataStateless<DisplacementDim> swelling_data;
    ElasticTangentStiffnessData<DisplacementDim> C_el_data;
    BiotData biot_data;
    SolidCompressibilityData solid_compressibility_data;
    SaturationDataDeriv dS_L_data;
    BishopsData bishops_data;
    // TODO why not usual state tracking for that?
    PrevState<BishopsData> bishops_data_prev;
    SolidThermalExpansionData<DisplacementDim> s_therm_exp_data;
    PermeabilityData<DisplacementDim> perm_data;
    FluidThermalExpansionData f_therm_exp_data;
    EquivalentPlasticStrainData equiv_plast_strain_data;
};
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
