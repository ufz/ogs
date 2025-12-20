// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
using StiffnessTensor = BaseLib::StrongType<KelvinMatrix<DisplacementDim>,
                                            struct StiffnessTensorTag>;
}  // namespace ProcessLib::RichardsMechanics
