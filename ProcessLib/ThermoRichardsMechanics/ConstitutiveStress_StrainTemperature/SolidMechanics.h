/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Biot.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Bishops.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EquivalentPlasticStrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/MaterialState.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidThermalExpansion.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Swelling.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TotalStressData.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
struct SolidMechanicsDataStateful
{
    // TODO it seems fragile that some data have to be initialized that way.
    KelvinVector<DisplacementDim> sigma_eff = KVzero<DisplacementDim>();
    KelvinVector<DisplacementDim> eps_m = KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = SolidMechanicsDataStateful<DisplacementDim>;

        // TODO add eps_m?
        return ProcessLib::Reflection::reflectWithName("sigma",
                                                       &Self::sigma_eff);
    }
};

template <int DisplacementDim>
struct SolidMechanicsDataStateless
{
    KelvinMatrix<DisplacementDim> stiffness_tensor = KMnan<DisplacementDim>();
    KelvinVector<DisplacementDim> J_uT_BT_K_N = KVnan<DisplacementDim>();
    KelvinVector<DisplacementDim> J_up_BT_K_N = KVnan<DisplacementDim>();
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
        PrevState<SolidMechanicsDataStateful<DisplacementDim>> const&
            prev_state,
        SolidMechanicsDataStateful<DisplacementDim>& current_state,
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
