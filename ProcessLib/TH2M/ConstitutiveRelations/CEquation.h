/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "Biot.h"
#include "ConstitutiveDensity.h"
#include "PhaseTransitionData.h"
#include "Porosity.h"
#include "Saturation.h"
#include "SolidCompressibility.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct FC2aData
{
    double a = nan;
};

struct FC2aDerivativeData
{
    double dp_GR = nan;
    double dp_cap = nan;
    double dT = nan;
};

struct FC2aModel
{
    void eval(BiotData const biot_data,
              CapillaryPressureData const pCap,
              ConstituentDensityData const& constituent_density_data,
              PorosityData const& porosity_data,
              SaturationData const& S_L_data,
              SolidCompressibilityData const beta_p_SR,
              FC2aData& fC_2a) const;

    void dEval(BiotData const& biot_data,
               CapillaryPressureData const pCap,
               ConstituentDensityData const& constituent_density_data,
               PhaseTransitionData const& phase_transition_data,
               PorosityData const& porosity_data,
               PorosityDerivativeData const& porosity_d_data,
               SaturationData const& S_L_data,
               SaturationDataDeriv const& dS_L_dp_cap,
               SolidCompressibilityData const& beta_p_SR,
               FC2aDerivativeData& dfC_2a) const;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
