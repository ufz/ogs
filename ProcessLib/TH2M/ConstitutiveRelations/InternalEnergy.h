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
#include "BaseLib/StrongType.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
using InternalEnergyData =
    BaseLib::StrongType<double, struct InternalEnergyTag>;

struct EffectiveVolumetricInternalEnergyDerivatives
{
    double drho_u_eff_dT = nan;
    double drho_u_eff_dp_GR = nan;
    double drho_u_eff_dp_cap = nan;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
