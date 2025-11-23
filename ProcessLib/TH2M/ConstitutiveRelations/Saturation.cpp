/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Saturation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
void SaturationModel::eval(SpaceTimeData const& x_t,
                           MediaData const& media_data,
                           CapillaryPressureData const& p_cap,
                           SaturationData& S_L_data) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.capillary_pressure = p_cap.pCap;

    S_L_data.S_L = media_data.saturation_prop.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);
}

void SaturationModel::dEval(SpaceTimeData const& x_t,
                            MediaData const& media_data,
                            CapillaryPressureData const& p_cap,
                            SaturationDataDeriv& dS_L_data) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.capillary_pressure = p_cap.pCap;

    dS_L_data() = media_data.saturation_prop.template dValue<double>(
        variables, MPL::Variable::capillary_pressure, x_t.x, x_t.t, x_t.dt);
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
