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
using ConstitutiveModels = std::tuple<
    ElasticTangentStiffnessModel<DisplacementDim>,
    BiotModel,
    SolidCompressibilityModel<DisplacementDim,
                              SolidConstitutiveRelation<DisplacementDim>>,
    SaturationModel<DisplacementDim>,
    BishopsModel,
    BishopsPrevModel,
    PorosityModel<DisplacementDim>,

    SwellingModel<DisplacementDim>,
    SolidThermalExpansionModel<DisplacementDim>,
    SolidMechanicsModel<DisplacementDim>,
    LiquidDensityModel<DisplacementDim>,

    SolidDensityModel<DisplacementDim>,
    GravityModel<DisplacementDim>,
    LiquidViscosityModel<DisplacementDim>,
    TransportPorosityModel<DisplacementDim>,
    PermeabilityModel<DisplacementDim>,
    ThermoOsmosisModel<DisplacementDim>,
    DarcyLawModel<DisplacementDim>,
    TRMHeatStorageAndFluxModel<DisplacementDim>,
    TRMVaporDiffusionModel<DisplacementDim>,
    FluidThermalExpansionModel<DisplacementDim>,
    TRMStorageModel<DisplacementDim>,
    EqPModel<DisplacementDim>,
    EqTModel<DisplacementDim>>;

template <int DisplacementDim, typename TRMProcessData>
ConstitutiveModels<DisplacementDim> createConstitutiveModels(
    TRMProcessData const& process_data,
    SolidConstitutiveRelation<DisplacementDim> const& solid_material)
{
    return {
        ElasticTangentStiffnessModel<DisplacementDim>{solid_material},
        {} /* BiotModel */,
        SolidCompressibilityModel<DisplacementDim,
                                  SolidConstitutiveRelation<DisplacementDim>>{
            solid_material},
        {} /* SaturationModel<DisplacementDim> */,
        {} /* BishopsModel */,
        {} /* BishopsPrevModel */,
        {} /* PorosityModel<DisplacementDim> */,

        {} /* SwellingModel<DisplacementDim> */,
        {} /* SolidThermalExpansionModel<DisplacementDim> */,
        SolidMechanicsModel<DisplacementDim>{solid_material},
        {} /* LiquidDensityModel<DisplacementDim> */,

        {} /* SolidDensityModel<DisplacementDim> */,
        GravityModel<DisplacementDim>{process_data.specific_body_force},
        {} /* LiquidViscosityModel<DisplacementDim> */,
        {} /* TransportPorosityModel<DisplacementDim> */,
        {} /* PermeabilityModel<DisplacementDim> */,
        {} /* ThermoOsmosisModel<DisplacementDim> */,
        DarcyLawModel<DisplacementDim>{process_data.specific_body_force},
        {} /* TRMHeatStorageAndFluxModel<DisplacementDim> */,
        {} /* TRMVaporDiffusionModel<DisplacementDim> */,
        {} /* FluidThermalExpansionModel<DisplacementDim> */,
        {} /* TRMStorageModel<DisplacementDim> */,
        EqPModel<DisplacementDim>{process_data.specific_body_force},
        {} /* EqTModel<DisplacementDim> */
    };
}
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
