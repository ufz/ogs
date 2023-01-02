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

#include "NumericalStabilization.h"

namespace NumLib
{
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

}  // namespace NumLib
