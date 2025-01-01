/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Base.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
using StiffnessTensor = BaseLib::StrongType<KelvinMatrix<DisplacementDim>,
                                            struct StiffnessTensorTag>;
}  // namespace ProcessLib::RichardsMechanics
