/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 24, 2022, 12:53 PM
 */

#pragma once

#include <algorithm>
#include <limits>
#include <numeric>
#include <range/v3/view/filter.hpp>
#include <typeinfo>

#include "NumericalStabilization.h"

namespace NumLib
{
template <typename IPData, typename FluxVectorType, typename Derived>
void assembleOriginalAdvectionMatrix(
    IPData const& ip_data_vector,
    FluxVectorType const& ip_flux_vector,
    Eigen::MatrixBase<Derived>& laplacian_matrix)
{
    for (std::size_t ip = 0; ip < ip_flux_vector.size(); ++ip)
    {
        auto const& ip_data = ip_data_vector[ip];
        auto const w = ip_data.integration_weight;
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        laplacian_matrix.noalias() +=
            N.transpose() * ip_flux_vector[ip].transpose() * dNdx * w;
    }
}

template <typename Derived>
void applyFullUpwind(Eigen::VectorXd const& quasi_nodal_flux,
                     Eigen::MatrixBase<Derived>& diffusion_matrix)
{
    std::vector<int> node_ids(quasi_nodal_flux.size());
    std::iota(std::begin(node_ids), std::end(node_ids), 0);
    auto down_wind = [&quasi_nodal_flux](const int id)
    { return quasi_nodal_flux[id] < 0; };
    auto up_wind = [&quasi_nodal_flux](const int id)
    { return quasi_nodal_flux[id] >= 0; };

    double q_in = 0.0;
    for (auto const downwind_node_id :
         node_ids | ranges::views::filter(down_wind))
    {
        q_in -= quasi_nodal_flux[downwind_node_id];
    }

    if (q_in < std::numeric_limits<double>::epsilon())
    {
        return;
    }

    for (auto const upwind_node_id : node_ids | ranges::views::filter(up_wind))
    {
        diffusion_matrix.diagonal()[upwind_node_id] +=
            quasi_nodal_flux[upwind_node_id];
    }

    for (auto const downwind_node_id :
         node_ids | ranges::views::filter(down_wind))
    {
        double const row_factor = quasi_nodal_flux[downwind_node_id] / q_in;
        for (auto const upwind_node_id :
             node_ids | ranges::views::filter(up_wind))
        {
            diffusion_matrix.row(downwind_node_id)[upwind_node_id] +=
                row_factor * quasi_nodal_flux[upwind_node_id];
        }
    }
}

template <typename IPData, typename FluxVectorType, typename Derived>
void applyFullUpwind(IPData const& ip_data_vector,
                     FluxVectorType const& ip_flux_vector,
                     Eigen::MatrixBase<Derived>& laplacian_matrix)
{
    Eigen::VectorXd quasi_nodal_flux(laplacian_matrix.rows());
    quasi_nodal_flux.setZero();

    for (std::size_t ip = 0; ip < ip_flux_vector.size(); ++ip)
    {
        auto const& ip_data = ip_data_vector[ip];
        auto const w = ip_data.integration_weight;
        auto const& dNdx = ip_data.dNdx;

        quasi_nodal_flux.noalias() -= ip_flux_vector[ip].transpose() * dNdx * w;
    }

    applyFullUpwind(quasi_nodal_flux, laplacian_matrix);
}

template <typename IPData, typename FluxVectorType, typename Derived>
void assembleAdvectionMatrix(NumericalStabilization const* const stabilizer,
                             IPData const& ip_data_vector,
                             double const average_velocity,
                             FluxVectorType const& ip_flux_vector,
                             Eigen::MatrixBase<Derived>& laplacian_matrix)
{
    if (stabilizer)
    {
        auto const& stabilizer_ref = *(stabilizer);
        if (typeid(stabilizer_ref) == typeid(NumLib::FullUpwind) &&
            (average_velocity > stabilizer->getCutoffVelocity()))
        {
            applyFullUpwind(ip_data_vector, ip_flux_vector, laplacian_matrix);
            return;
        }
    }

    assembleOriginalAdvectionMatrix(ip_data_vector, ip_flux_vector,
                                    laplacian_matrix);
}

}  // namespace NumLib
