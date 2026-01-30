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
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct TemperatureData<2>;
extern template struct TemperatureData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
