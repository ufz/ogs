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

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct SolidMechanicsDataStateless
{
    KelvinMatrix<DisplacementDim> stiffness_tensor = KMnan<DisplacementDim>();
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
