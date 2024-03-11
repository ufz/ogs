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

#include "ElasticTangentStiffnessModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/DarcyLaw.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqP.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqT.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/FluidThermalExpansion.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Gravity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/PermeabilityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidDensity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMHeatStorageAndFlux.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMStorage.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMVaporDiffusion.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStress_StrainTemperature/SolidCompressibilityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStress_StrainTemperature/SolidMechanics.h"
#include "Swelling.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
/// Constitutive models used for assembly.
template <int DisplacementDim>
struct ConstitutiveModels
{
    template <typename TRMProcessData>
    explicit ConstitutiveModels(
        TRMProcessData const& process_data,
        SolidConstitutiveRelation<DisplacementDim> const& solid_material)
        : elastic_tangent_stiffness_model(solid_material),
          solid_compressibility_model(solid_material),
          s_mech_model(solid_material),
          grav_model(process_data.specific_body_force),
          darcy_model(process_data.specific_body_force),
          eq_p_model(process_data.specific_body_force)
    {
    }

    ElasticTangentStiffnessModel<DisplacementDim>
        elastic_tangent_stiffness_model;
    BiotModel biot_model;
    SolidCompressibilityModel<DisplacementDim,
                              SolidConstitutiveRelation<DisplacementDim>>
        solid_compressibility_model;
    SaturationModel<DisplacementDim> S_L_model;
    BishopsModel bishops_model;
    BishopsPrevModel bishops_prev_model;
    PorosityModel<DisplacementDim> poro_model;
    TransportPorosityModel<DisplacementDim> transport_poro_model;
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
    EqPModel<DisplacementDim> eq_p_model;
    EqTModel<DisplacementDim> eq_T_model;
    ThermoOsmosisModel<DisplacementDim> th_osmosis_model;
};
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
