// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/Graph/ConstructModels.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/BishopsModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/DarcyLawModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqPModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EqTModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/FluidThermalExpansionModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/GravityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/LiquidDensityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/LiquidViscosityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/PermeabilityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/PorosityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidDensityModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidThermalExpansion.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMHeatStorageAndFluxModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMStorageModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TRMVaporDiffusionModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/ThermoOsmosisModel.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStressSaturation_StrainPressureTemperature/EffectiveStressModel.h"
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
    EffectiveStressModel<DisplacementDim>,
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
        ProcessLib::ConstitutiveRelations::SpecificBodyForce<DisplacementDim>(
            process_data.specific_body_force),
        solid_material);
}

}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
