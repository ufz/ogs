/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "AverageMolarMass.h"
#include "BishopsPowerLaw.h"
#include "BishopsSaturationCutoff.h"
#include "CapillaryPressureSaturation/SaturationBrooksCorey.h"
#include "CapillaryPressureSaturation/SaturationLiakopoulos.h"
#include "CapillaryPressureSaturation/SaturationVanGenuchten.h"
#include "ClausiusClapeyron.h"
#include "Constant.h"
#include "Curve.h"
#include "Density/WaterVapourDensity.h"
#include "DupuitPermeability.h"
#include "EmbeddedFracturePermeability.h"
#include "Enthalpy/LinearWaterVapourLatentHeat.h"
#include "Exponential.h"
#include "GasPressureDependentPermeability.h"
#include "IdealGasLaw.h"
#include "Linear.h"
#include "Parameter.h"
#include "PorosityFromMassBalance.h"
#include "RelativePermeability/RelPermBrooksCorey.h"
#include "RelativePermeability/RelPermBrooksCoreyNonwettingPhase.h"
#include "RelativePermeability/RelPermLiakopoulos.h"
#include "RelativePermeability/RelPermNonWettingPhaseVanGenuchtenMualem.h"
#include "RelativePermeability/RelPermUdell.h"
#include "RelativePermeability/RelPermUdellNonwettingPhase.h"
#include "RelativePermeability/RelPermVanGenuchten.h"
#include "SaturationDependentHeatConduction.h"
#include "SaturationDependentSwelling.h"
#include "ThermalConductivity/SoilThermalConductivitySomerton.h"
#include "TransportPorosityFromMassBalance.h"
#include "VapourDiffusion/VapourDiffusionFEBEX.h"
#include "VapourDiffusion/VapourDiffusionPMQ.h"
