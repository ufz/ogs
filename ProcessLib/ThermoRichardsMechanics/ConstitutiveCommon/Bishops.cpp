// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Bishops.h"

namespace ProcessLib::ThermoRichardsMechanics
{
static void bishopsModelEvalImpl(SpaceTimeData const& x_t,
                                 MediaData const& media_data,
                                 SaturationData const& S_L_data,
                                 BishopsData& out)
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.liquid_saturation = S_L_data.S_L;

    auto const& medium = media_data.medium;

    out.chi_S_L = medium.property(MPL::PropertyType::bishops_effective_stress)
                      .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    out.dchi_dS_L = medium.property(MPL::PropertyType::bishops_effective_stress)
                        .template dValue<double>(
                            variables, MPL::Variable::liquid_saturation, x_t.x,
                            x_t.t, x_t.dt);
}

void BishopsModel::eval(SpaceTimeData const& x_t, MediaData const& media_data,
                        SaturationData const& S_L_data, BishopsData& out) const
{
    bishopsModelEvalImpl(x_t, media_data, S_L_data, out);
}

void BishopsPrevModel::eval(SpaceTimeData const& x_t,
                            MediaData const& media_data,
                            PrevState<SaturationData> const& S_L_data,
                            PrevState<BishopsData>& out) const
{
    bishopsModelEvalImpl(x_t, media_data, *S_L_data, *out);
}
}  // namespace ProcessLib::ThermoRichardsMechanics
