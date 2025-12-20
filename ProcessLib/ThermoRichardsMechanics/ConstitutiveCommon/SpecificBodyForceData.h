// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SpecificBodyForceData
{
    Eigen::Vector<double, DisplacementDim> specific_body_force;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
