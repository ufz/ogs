/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BaseLib/StrongType.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
using EquivalentPlasticStrainData =
    BaseLib::StrongType<double, struct EquivalentPlasticStrainTag>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
