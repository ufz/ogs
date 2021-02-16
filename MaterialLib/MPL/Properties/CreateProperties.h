/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "CapillaryPressureSaturation/CreateCapillaryPressureRegularizedVanGenuchten.h"
#include "CapillaryPressureSaturation/CreateCapillaryPressureVanGenuchten.h"
#include "CapillaryPressureSaturation/CreateSaturationBrooksCorey.h"
#include "CapillaryPressureSaturation/CreateSaturationLiakopoulos.h"
#include "CapillaryPressureSaturation/CreateSaturationVanGenuchten.h"
#include "CreateAverageMolarMass.h"
#include "CreateBishopsPowerLaw.h"
#include "CreateBishopsSaturationCutoff.h"
#include "CreateClausiusClapeyron.h"
#include "CreateConstant.h"
#include "CreateCurve.h"
#include "CreateDupuitPermeability.h"
#include "CreateEmbeddedFracturePermeability.h"
#include "CreateExponential.h"
#include "CreateIdealGasLaw.h"
#include "CreateKozenyCarmanModel.h"
#include "CreateLinear.h"
#include "CreateParameter.h"
#include "CreatePermeabilityMohrCoulombFailureIndexModel.h"
#include "CreatePermeabilityOrthotropicPowerLaw.h"
#include "CreatePorosityFromMassBalance.h"
#include "CreateSaturationDependentHeatConduction.h"
#include "CreateSaturationDependentSwelling.h"
#include "CreateStrainDependentPermeability.h"
#include "CreateTransportPorosityFromMassBalance.h"
#include "RelativePermeability/CreateRelPermBrooksCorey.h"
#include "RelativePermeability/CreateRelPermBrooksCoreyNonwettingPhase.h"
#include "RelativePermeability/CreateRelPermLiakopoulos.h"
#include "RelativePermeability/CreateRelPermNonWettingPhaseVanGenuchtenMualem.h"
#include "RelativePermeability/CreateRelPermUdell.h"
#include "RelativePermeability/CreateRelPermUdellNonwettingPhase.h"
#include "RelativePermeability/CreateRelPermVanGenuchten.h"
#include "SwellingStress/CreateLinearSaturationSwellingStress.h"
