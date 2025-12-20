// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <vector>

#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"

namespace NumLib
{
/**
 * The function computes the element average of \f$\nabla N_i\f$ by
 * \f[
 *     \bar{\nabla N_i} = \frac{\int_{\Omega_e}\nabla N_i{\mathrm d}\Omega}
 *                        {\int_{\Omega_e}{\mathrm d}\Omega}.
 * \f]
 * @param local_node_id        Node ID in an element,
 * @param element              Element,
 * @param integration_method   Integration method,
 * @param ip_data              Integration point data,
 * @param is_axially_symmetric Indicator for axisymmetry.
 * @return Averaged \f$\nabla N_i\f$ of an element.
 */
template <int DisplacementDim,
          typename ShapeFunction,
          typename ShapeMatricesType,
          typename IpData>
Eigen::Vector3d averageGradShapeFunction(
    int const local_node_id,
    MeshLib::Element const& element,
    NumLib::GenericIntegrationMethod const& integration_method,
    std::vector<IpData, Eigen::aligned_allocator<IpData>> const& ip_data,
    const bool is_axially_symmetric)
{
    Eigen::Vector3d bar_gradN = Eigen::Vector3d::Zero();
    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();
    assert(n_integration_points == ip_data.size());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = ip_data[ip].N_u;
        auto const& dNdx = ip_data[ip].dNdx_u;

        auto const dNidx = dNdx.col(local_node_id);

        auto const& w = ip_data[ip].integration_weight;
        bar_gradN.template segment<DisplacementDim>(0) += dNidx * w;

        if (is_axially_symmetric)
        {
            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(element, N);
            bar_gradN[2] += w * N(local_node_id) / x_coord;
        }
    }
    return bar_gradN;
}
}  // namespace NumLib
