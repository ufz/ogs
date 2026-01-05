// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
template <typename ShapeMatricesType>
struct IntegrationPointData final
{
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    typename ShapeMatricesType::GlobalDimVectorType v_darcy;

    double saturation = std::numeric_limits<double>::quiet_NaN();
    double saturation_prev = std::numeric_limits<double>::quiet_NaN();
    double porosity = std::numeric_limits<double>::quiet_NaN();
    double porosity_prev = std::numeric_limits<double>::quiet_NaN();
    double transport_porosity = std::numeric_limits<double>::quiet_NaN();
    double transport_porosity_prev = std::numeric_limits<double>::quiet_NaN();
    double dry_density_solid = std::numeric_limits<double>::quiet_NaN();
    double dry_density_pellet_saturated =
        std::numeric_limits<double>::quiet_NaN();
    double dry_density_pellet_unsaturated =
        std::numeric_limits<double>::quiet_NaN();

    double integration_weight;

    void pushBackState()
    {
        saturation_prev = saturation;
        porosity_prev = porosity;
        transport_porosity_prev = transport_porosity;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
