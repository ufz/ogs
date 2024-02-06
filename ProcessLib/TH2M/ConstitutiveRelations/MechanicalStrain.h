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
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "SolidThermalExpansion.h"
#include "Swelling.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct MechanicalStrainData
{
    // TODO it seems fragile that some data have to be initialized that way.
    KelvinVector<DisplacementDim> eps_m = KV::KVzero<DisplacementDim>();
};

template <int DisplacementDim>
struct MechanicalStrainModel
{
    void eval(
        TemperatureData const& T_data,
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        ProcessLib::ConstitutiveRelations::StrainData<DisplacementDim> const&
            strain_data,
        KelvinVector<DisplacementDim> const& eps_prev,
        PrevState<MechanicalStrainData<DisplacementDim>> const& eps_m_prev,
        SwellingDataStateless<DisplacementDim> const& swelling_data,
        MechanicalStrainData<DisplacementDim>& out) const;
};

extern template struct MechanicalStrainModel<2>;
extern template struct MechanicalStrainModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
