/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>

#include "Constitutive/DarcyLaw.h"
#include "Constitutive/EqP.h"
#include "Constitutive/EqT.h"
#include "Constitutive/EqU.h"
#include "Constitutive/FluidThermalExpansion.h"
#include "Constitutive/Gravity.h"
#include "Constitutive/MaterialState.h"
#include "Constitutive/Porosity.h"
#include "Constitutive/Saturation.h"
#include "Constitutive/SolidMechanics.h"
#include "Constitutive/TRMHeatStorageAndFlux.h"
#include "Constitutive/TRMStorage.h"
#include "Constitutive/TRMVaporDiffusion.h"
#include "Constitutive/ThermoOsmosis.h"
#include "ThermoRichardsMechanicsProcessData.h"

namespace ProcessLib::ThermoRichardsMechanics
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

        return ProcessLib::Reflection::reflectWithoutName(
            &Self::S_L_data,
            &Self::poro_data,
            &Self::transport_poro_data,
            &Self::eps_data,
            &Self::swelling_data,
            &Self::s_mech_data);
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

        return ProcessLib::Reflection::reflectWithoutName(&Self::darcy_data,
                                                          &Self::rho_L_data,
                                                          &Self::mu_L_data,
                                                          &Self::rho_S_data);
    }
};

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
struct ConstitutiveData
{
    SwellingDataStateless<DisplacementDim> swelling_data;
    SolidMechanicsDataStateless<DisplacementDim> s_mech_data;
    GravityData<DisplacementDim> grav_data;
    TRMHeatStorageAndFluxData<DisplacementDim> heat_data;
    TRMVaporDiffusionData<DisplacementDim> vap_data;
    TRMStorageData storage_data;
    EqUData<DisplacementDim> eq_u_data;
    EqPData<DisplacementDim> eq_p_data;
    EqTData<DisplacementDim> eq_T_data;
    ThermoOsmosisData<DisplacementDim> th_osmosis_data;
};

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
struct ConstitutiveTempData
{
    ElasticTangentStiffnessData<DisplacementDim> C_el_data;
    BiotData biot_data;
    SolidCompressibilityData solid_compressibility_data;
    SaturationDataDeriv dS_L_data;
    BishopsData bishops_data;
    // TODO why not usual state tracking for that?
    BishopsData bishops_data_prev;
    SolidThermalExpansionData<DisplacementDim> s_therm_exp_data;
    PermeabilityData<DisplacementDim> perm_data;
    FluidThermalExpansionData f_therm_exp_data;
};

/// Constitutive models used for assembly.
template <int DisplacementDim>
struct ConstitutiveModels
{
    explicit ConstitutiveModels(
        ThermoRichardsMechanicsProcessData<DisplacementDim> const& process_data,
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : elastic_tangent_stiffness_model(solid_material),
          solid_compressibility_model(solid_material),
          s_mech_model(solid_material),
          grav_model(process_data.specific_body_force),
          darcy_model(process_data.specific_body_force),
          eq_u_model(process_data.specific_body_force),
          eq_p_model(process_data.specific_body_force)
    {
    }

    ElasticTangentStiffnessModel<DisplacementDim>
        elastic_tangent_stiffness_model;
    BiotModel biot_model;
    SolidCompressibilityModel<DisplacementDim> solid_compressibility_model;
    SaturationModel<DisplacementDim> S_L_model;
    BishopsModel bishops_model;
    PorosityModel<DisplacementDim> poro_model;
    SwellingModel<DisplacementDim> swelling_model;
    SolidThermalExpansionModel<DisplacementDim> s_therm_exp_model;
    SolidMechanicsModel<DisplacementDim> s_mech_model;
    LiquidDensityModel<DisplacementDim> rho_L_model;
    SolidDensityModel<DisplacementDim> rho_S_model;
    GravityModel<DisplacementDim> grav_model;
    LiquidViscosityModel<DisplacementDim> mu_L_model;
    PermeabilityModel<DisplacementDim> perm_model;
    DarcyLawModel<DisplacementDim> darcy_model;
    TRMHeatStorageAndFluxModel<DisplacementDim> heat_storage_and_flux_model;
    TRMVaporDiffusionModel<DisplacementDim> vapor_diffusion_model;
    FluidThermalExpansionModel<DisplacementDim> f_therm_exp_model;
    TRMStorageModel<DisplacementDim> storage_model;
    EqUModel<DisplacementDim> eq_u_model;
    EqPModel<DisplacementDim> eq_p_model;
    EqTModel<DisplacementDim> eq_T_model;
    ThermoOsmosisModel<DisplacementDim> th_osmosis_model;
};

template <int DisplacementDim>
struct ConstitutiveSetting
{
    explicit ConstitutiveSetting(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        ThermoRichardsMechanicsProcessData<DisplacementDim> const& process_data)
        : process_data_(process_data), solid_material_(solid_material)
    {
    }

    /// Evaluate the constitutive setting.
    void eval(ConstitutiveModels<DisplacementDim>& models, double const t,
              double const dt, ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium& medium,
              TemperatureData<DisplacementDim> const& T_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              KelvinVector<DisplacementDim> const& eps_arg,
              KelvinVector<DisplacementDim> const& eps_prev_arg,
              StatefulData<DisplacementDim>& state,
              StatefulData<DisplacementDim> const& prev_state,
              MaterialStateData<DisplacementDim>& mat_state,
              ConstitutiveTempData<DisplacementDim>& tmp,
              OutputData<DisplacementDim>& out,
              ConstitutiveData<DisplacementDim>& cd);

private:
    ThermoRichardsMechanicsProcessData<DisplacementDim> const& process_data_;
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
