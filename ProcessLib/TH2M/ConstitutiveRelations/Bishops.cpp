// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Bishops.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
static void bishopsModelEvalImpl(SpaceTimeData const& x_t,
                                 MediaData const& media_data,
                                 SaturationData const& S_L_data,
                                 BishopsData& out)
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.liquid_saturation = S_L_data.S_L;

    out.chi_S_L =
        media_data.bishops_effective_stress_prop.template value<double>(
            variables, x_t.x, x_t.t, x_t.dt);

    out.dchi_dS_L =
        media_data.bishops_effective_stress_prop.template dValue<double>(
            variables, MPL::Variable::liquid_saturation, x_t.x, x_t.t, x_t.dt);
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

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
