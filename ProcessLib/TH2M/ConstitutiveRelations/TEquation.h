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
#include "InternalEnergy.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct FT1Data
{
    double m = nan;
};

struct FT1DerivativeData
{
    double dp_GR = nan;
    double dp_cap = nan;
    double dT = nan;
};

struct FT1Model
{
    void eval(double const dt,
              InternalEnergyData const& internal_energy_data,
              PrevState<InternalEnergyData> const& internal_energy_data_prev,
              FT1Data& fT_1) const;

    void dEval(double const dt,
               EffectiveVolumetricInternalEnergyDerivatives const&
                   effective_volumetric_internal_energy_d_data,
               FT1DerivativeData& dfT_1) const;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
