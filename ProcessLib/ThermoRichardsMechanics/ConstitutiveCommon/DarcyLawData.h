// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <string_view>

#include "BaseLib/StrongType.h"
namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
using DarcyLawData = BaseLib::StrongType<Eigen::Vector<double, DisplacementDim>,
                                         struct DarcyLawDataTag>;

constexpr std::string_view ioName(struct DarcyLawDataTag*)
{
    return "velocity";
}

}  // namespace ProcessLib::ThermoRichardsMechanics