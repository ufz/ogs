// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ProcessLib/ConstitutiveRelations/MechanicalStrainData.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "SolidThermalExpansion.h"
#include "Swelling.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

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
