/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ConstitutiveRelations/Base.h"

namespace ProcessLib::RichardsMechanics
{
// TODO directly declare these type aliases in Traits.h
/// Data whose state must be tracked by the TRM process.
template <int DisplacementDim>
using StatefulData = std::tuple<>;

template <int DisplacementDim>
using StatefulDataPrev = ProcessLib::ConstitutiveRelations::PrevStateOf<
    StatefulData<DisplacementDim>>;

/// Data that is needed for output purposes, but not directly for the assembly.
template <int DisplacementDim>
using OutputData = std::tuple<>;

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
using ConstitutiveData = std::tuple<>;

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
using ConstitutiveTempData = std::tuple<>;
}  // namespace ProcessLib::RichardsMechanics
