// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Biot.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
void BiotModel::eval(SpaceTimeData const& x_t, MediaData const& media_data,
                     BiotData& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;

    *out = media_data.biot_coefficient_prop.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
