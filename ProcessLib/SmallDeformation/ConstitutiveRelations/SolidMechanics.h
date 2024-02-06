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
#include "FreeEnergyDensity.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialState.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ConstitutiveRelations/StressData.h"

namespace ProcessLib::SmallDeformation
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

    void eval(SpaceTimeData const& x_t,
              Temperature const& temperature,
              StrainData<DisplacementDim> const& eps_data,
              PrevState<StrainData<DisplacementDim>> const& eps_data_prev,
              MaterialStateData<DisplacementDim>& mat_state,
              PrevState<StressData<DisplacementDim>> const& stress_data_prev,
              StressData<DisplacementDim>& stress_data,
              SolidMechanicsDataStateless<DisplacementDim>& current_stateless,
              FreeEnergyDensityData& free_energy_density_data) const;

    auto getInternalVariables() const
    {
        return solid_material_.getInternalVariables();
    }

private:
    SolidConstitutiveRelation<DisplacementDim> const& solid_material_;
};

extern template struct SolidMechanicsModel<2>;
extern template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
