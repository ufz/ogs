/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
/// Real constituent partial densities.
struct ConstituentDensityData
{
    // Gas phase density
    double rho_C_GR = nan;
    double rho_W_GR = nan;

    // Liquid phase density
    double rho_C_LR = nan;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
