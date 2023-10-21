/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "AverageMolarMass.h"
#include "BishopsPowerLaw.h"
#include "BishopsSaturationCutoff.h"
#include "CapillaryPressureSaturation/SaturationBrooksCorey.h"
#include "CapillaryPressureSaturation/SaturationExponential.h"
#include "CapillaryPressureSaturation/SaturationLiakopoulos.h"
#include "CapillaryPressureSaturation/SaturationVanGenuchten.h"
#include "CapillaryPressureSaturation/SaturationVanGenuchtenWithVolumetricStrain.h"
#include "ClausiusClapeyron.h"
#include "Constant.h"
#include "Curve.h"
#include "Density/WaterDensityIAPWSIF97Region1.h"
#include "Density/WaterVapourDensity.h"
#include "DupuitPermeability.h"
#include "EffectiveThermalConductivityPorosityMixing.h"
#include "EmbeddedFracturePermeability.h"
#include "Enthalpy/LinearWaterVapourLatentHeat.h"
#include "Enthalpy/WaterVapourLatentHeatWithCriticalTemperature.h"
#include "Exponential.h"
#include "Function.h"
#include "GasPressureDependentPermeability.h"
#include "IdealGasLaw.h"
#include "IdealGasLawBinaryMixture.h"
#include "Linear.h"
#include "OrthotropicEmbeddedFracturePermeability.h"
#include "Parameter.h"
#include "PorosityFromMassBalance.h"
#include "RelativePermeability/RelPermBrooksCorey.h"
#include "RelativePermeability/RelPermBrooksCoreyNonwettingPhase.h"
#include "RelativePermeability/RelPermGeneralizedPower.h"
#include "RelativePermeability/RelPermGeneralizedPowerNonwettingPhase.h"
#include "RelativePermeability/RelPermLiakopoulos.h"
#include "RelativePermeability/RelPermNonWettingPhaseVanGenuchtenMualem.h"
#include "RelativePermeability/RelPermUdell.h"
#include "RelativePermeability/RelPermUdellNonwettingPhase.h"
#include "RelativePermeability/RelPermVanGenuchten.h"
#include "SaturationDependentSwelling.h"
#include "SpecificHeatCapacityWithLatentHeat.h"
#include "TemperatureDependentDiffusion.h"
#include "TemperatureDependentFraction.h"
#include "ThermalConductivity/SaturationWeightedThermalConductivity.h"
#include "TransportPorosityFromMassBalance.h"
#include "VapourDiffusion/VapourDiffusionDeVries.h"
#include "VapourDiffusion/VapourDiffusionFEBEX.h"
#include "VapourDiffusion/VapourDiffusionPMQ.h"
#include "Viscosity/LiquidViscosityVogels.h"
#include "Viscosity/WaterViscosityIAPWS.h"
#include "VolumeFractionAverage.h"
