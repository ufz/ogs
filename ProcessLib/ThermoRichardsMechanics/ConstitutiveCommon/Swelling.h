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

#include "ElasticTangentStiffnessData.h"
#include "Saturation.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SwellingDataStateful
{
    KelvinVector<DisplacementDim> sigma_sw = KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = SwellingDataStateful<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("swelling_stress",
                                                       &Self::sigma_sw);
    }
};

template <int DisplacementDim>
struct SwellingDataStateless
{
    // TODO find a better name. Maybe swelling strain?
    KelvinVector<DisplacementDim> eps_m;
    KelvinVector<DisplacementDim> J_up_BT_K_N;
};

template <int DisplacementDim>
struct SwellingModel
{
    void eval(
        SpaceTimeData const& x_t, MediaData const& media_data,
        ElasticTangentStiffnessData<DisplacementDim> const& C_el_data,
        StrainData<DisplacementDim> const& eps_data,
        PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
        SaturationData const& S_L_data, SaturationDataDeriv const& dS_L_data,
        PrevState<SaturationData> const& S_L_prev_data,
        PrevState<SwellingDataStateful<DisplacementDim>> const& prev_state,
        SwellingDataStateful<DisplacementDim>& state,
        SwellingDataStateless<DisplacementDim>& out) const;
};

extern template struct SwellingModel<2>;
extern template struct SwellingModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
