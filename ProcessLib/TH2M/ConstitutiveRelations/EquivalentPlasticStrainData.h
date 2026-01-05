// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
