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
#include "EquivalentPlasticStrainData.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialState.h"
#include "MechanicalStrain.h"
#include "ProcessLib/ConstitutiveRelations/StressData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct SolidMechanicsDataStateless
{
    KelvinMatrix<DisplacementDim> stiffness_tensor =
        KV::KMnan<DisplacementDim>();
};

template <int DisplacementDim>
using SolidConstitutiveRelation =
    MaterialLib::Solids::MechanicsBase<DisplacementDim>;

template <int DisplacementDim>
struct SolidMechanicsModel
{
    explicit SolidMechanicsModel(
        SolidConstitutiveRelation<DisplacementDim> const& solid_material)
        : solid_material_(solid_material)
    {
    }

    void eval(
        const SpaceTimeData& x_t,
        TemperatureData const& T_data,
        MechanicalStrainData<DisplacementDim> const& mechanical_strain_data,
        PrevState<MechanicalStrainData<DisplacementDim>> const&
            mechanical_strain_prev_data,
        PrevState<ProcessLib::ConstitutiveRelations::StressData<
            DisplacementDim>> const& eff_stress_prev_data,
        ProcessLib::ConstitutiveRelations::StressData<DisplacementDim>&
            eff_stress_data,
        MaterialStateData<DisplacementDim>& mat_state,
        SolidMechanicsDataStateless<DisplacementDim>& out,
        EquivalentPlasticStrainData& equivalent_plastic_strain) const;

    auto getInternalVariables() const
    {
        return solid_material_.getInternalVariables();
    }

private:
    SolidConstitutiveRelation<DisplacementDim> const& solid_material_;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
