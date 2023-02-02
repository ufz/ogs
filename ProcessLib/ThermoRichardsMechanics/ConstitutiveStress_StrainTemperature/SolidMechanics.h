/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MechanicalStrainData.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Biot.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Bishops.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EquivalentPlasticStrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/MaterialState.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidMechanicsDataStateless.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidThermalExpansion.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TotalStressData.h"
#include "Swelling.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
struct EffectiveStressData
{
    // TODO it seems fragile that some data have to be initialized that way.
    KelvinVector<DisplacementDim> sigma_eff = KV::KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = EffectiveStressData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("sigma",
                                                       &Self::sigma_eff);
    }
};

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
        SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
        SwellingDataStateless<DisplacementDim> const& swelling_data,
        TemperatureData<DisplacementDim> const& T_data,
        CapillaryPressureData<DisplacementDim> const& p_cap_data,
        BiotData const& biot_data,
        BishopsData const& bishops_data,
        SaturationDataDeriv const& dS_L_data,
        StrainData<DisplacementDim> const& eps_data,
        PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
        MaterialStateData<DisplacementDim>& mat_state,
        PrevState<EffectiveStressData<DisplacementDim>> const&
            sigma_eff_prev_data,
        EffectiveStressData<DisplacementDim>& sigma_eff_data,
        PrevState<MechanicalStrainData<DisplacementDim>> const& eps_m_prev_data,
        MechanicalStrainData<DisplacementDim>& eps_m_data,
        TotalStressData<DisplacementDim>& total_stress_data,
        EquivalentPlasticStrainData& equiv_plast_strain_data,
        SolidMechanicsDataStateless<DisplacementDim>& out) const;

    auto getInternalVariables() const
    {
        return solid_material_.getInternalVariables();
    }

private:
    SolidConstitutiveRelation<DisplacementDim> const& solid_material_;
};

extern template struct SolidMechanicsModel<2>;
extern template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
