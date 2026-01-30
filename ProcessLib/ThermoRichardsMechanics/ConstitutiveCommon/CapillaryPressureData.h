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

// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct CapillaryPressureData<2>;
extern template struct CapillaryPressureData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
