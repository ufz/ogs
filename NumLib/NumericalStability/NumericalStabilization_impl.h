/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 24, 2022, 12:53 PM
 */

#pragma once

#include <limits>

#include "NumericalStabilization.h"

namespace NumLib
{
template <typename Derived>
void applyFullUpwind(Eigen::VectorXd const& quasi_nodal_flux,
                     Eigen::MatrixBase<Derived>& diffusion_matrix)
{
    double q_in = 0.0;
    std::vector<int> downwind_node_ids;
    std::vector<int> upwind_node_ids;
    for (int i = 0; i < quasi_nodal_flux.size(); i++)
    {
        double const q_i = quasi_nodal_flux[i];
        if (q_i < 0.0)
        {
            q_in -= q_i;
            downwind_node_ids.push_back(i);
        }
        else
        {
            upwind_node_ids.push_back(i);
        }
    }

    if (q_in < std::numeric_limits<double>::epsilon())
    {
        return;
    }

    for (int const upwind_node_id : upwind_node_ids)
    {
        diffusion_matrix.diagonal()[upwind_node_id] +=
            quasi_nodal_flux[upwind_node_id];
    }

    for (int const downwind_node_id : downwind_node_ids)
    {
        double const row_factor = quasi_nodal_flux[downwind_node_id] / q_in;
        for (int const upwind_node_id : upwind_node_ids)
        {
            diffusion_matrix.row(downwind_node_id)[upwind_node_id] +=
                row_factor * quasi_nodal_flux[upwind_node_id];
        }
    }
}

}  // namespace NumLib
