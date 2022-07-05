/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Biot.h"
#include "Bishops.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialState.h"
#include "SolidThermalExpansion.h"
#include "Swelling.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SolidMechanicsDataStateful
{
    // TODO it seems fragile that some data have to be initialized that way.
    KelvinVector<DisplacementDim> sigma_eff = KVzero<DisplacementDim>();
    KelvinVector<DisplacementDim> eps_m = KVzero<DisplacementDim>();
};

template <int DisplacementDim>
struct SolidMechanicsDataStateless
{
    KelvinVector<DisplacementDim> sigma_total = KVnan<DisplacementDim>();
    KelvinMatrix<DisplacementDim> stiffness_tensor = KMnan<DisplacementDim>();
    KelvinVector<DisplacementDim> J_uT_BT_K_N = KVnan<DisplacementDim>();
    double equivalent_plastic_strain = nan;
};

template <int DisplacementDim>
struct SolidMechanicsModel
{
    explicit SolidMechanicsModel(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
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
        StrainData<DisplacementDim> const& eps_data,
        StrainData<DisplacementDim> const& eps_prev_data,
        MaterialStateData<DisplacementDim>& mat_state,
        SolidMechanicsDataStateful<DisplacementDim> const& prev_state,
        SolidMechanicsDataStateful<DisplacementDim>& current_state,
        SolidMechanicsDataStateless<DisplacementDim>& out) const;

private:
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;
};

extern template struct SolidMechanicsModel<2>;
extern template struct SolidMechanicsModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
