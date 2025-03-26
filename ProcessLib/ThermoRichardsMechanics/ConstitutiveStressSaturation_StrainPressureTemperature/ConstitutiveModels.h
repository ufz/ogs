/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Graph/ConstructModels.h"
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
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStressSaturation_StrainPressureTemperature/SolidCompressibilityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStressSaturation_StrainPressureTemperature/SolidMechanics.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStressSaturation_StrainPressureTemperature/TransportPorosity.h"


namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
/// Constitutive models used for assembly.
template <int DisplacementDim>
using ConstitutiveModels = std::tuple<
    BiotModel,
    SolidMechanicsModel<DisplacementDim>,
    SolidCompressibilityModel<DisplacementDim,
                              SolidConstitutiveRelation<DisplacementDim>>,
    BishopsModel,
    BishopsPrevModel,
    PorosityModel<DisplacementDim>,

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

    SolidThermalExpansionModel<DisplacementDim>,
    FluidThermalExpansionModel<DisplacementDim>,
    TRMStorageModel<DisplacementDim>,
    EqPModel<DisplacementDim>,
    EqTModel<DisplacementDim>>;

template <int DisplacementDim, typename TRMProcessData>
ConstitutiveModels<DisplacementDim> createConstitutiveModels(
    TRMProcessData const& process_data,
    SolidConstitutiveRelation<DisplacementDim> const& solid_material)
{
    return ProcessLib::Graph::constructModels<
        ConstitutiveModels<DisplacementDim>>(
        SpecificBodyForceData<DisplacementDim>{
            process_data.specific_body_force},
        solid_material);
}

}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
