// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct CapillaryPressureData
{
    double p_cap;
    double p_cap_prev;
    Eigen::Vector<double, DisplacementDim> grad_p_cap;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
