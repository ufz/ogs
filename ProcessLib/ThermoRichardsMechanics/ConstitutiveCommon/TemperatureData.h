// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Dense>

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct TemperatureData
{
    double T;
    double T_prev;
    Eigen::Vector<double, DisplacementDim> grad_T;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
