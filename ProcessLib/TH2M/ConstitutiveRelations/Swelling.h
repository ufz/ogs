/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "ElasticTangentStiffnessData.h"
#include "Saturation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct SwellingDataStateful
{
    KelvinVector<DisplacementDim> sigma_sw = KV::KVzero<DisplacementDim>();

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
};

template <int DisplacementDim>
struct SwellingModel
{
    void eval(
        SpaceTimeData const& x_t, MediaData const& media_data,
        ElasticTangentStiffnessData<DisplacementDim> const& C_el_data,
        SaturationData const& S_L_data,
        PrevState<SaturationData> const& S_L_prev_data,
        PrevState<SwellingDataStateful<DisplacementDim>> const& prev_state,
        SwellingDataStateful<DisplacementDim>& state,
        SwellingDataStateless<DisplacementDim>& out) const;
};

extern template struct SwellingModel<2>;
extern template struct SwellingModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
