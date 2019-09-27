/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 7, 2016, 12:14 PM
 */

#include "CoupledSolutionsForStaggeredScheme.h"

#include "MathLib/LinAlg/LinAlg.h"
#include "Process.h"

namespace ProcessLib
{
CoupledSolutionsForStaggeredScheme::CoupledSolutionsForStaggeredScheme(
    std::vector<GlobalVector*> const& coupled_xs_)
    : coupled_xs(coupled_xs_)
{
    for (auto const* coupled_x : coupled_xs)
    {
        MathLib::LinAlg::setLocalAccessibleVector(*coupled_x);
    }
}

std::vector<std::vector<double>> getPreviousLocalSolutions(
    const CoupledSolutionsForStaggeredScheme& cpl_xs,
    const std::vector<std::vector<GlobalIndexType>>& indices)
{
    if (cpl_xs.coupled_xs_t0.empty())
    {
        return {};
    }

    const auto number_of_coupled_solutions = cpl_xs.coupled_xs.size();
    std::vector<std::vector<double>> local_xs_t0;
    local_xs_t0.reserve(number_of_coupled_solutions);

    int coupling_id = 0;
    for (auto const& x_t0 : cpl_xs.coupled_xs_t0)
    {
        local_xs_t0.emplace_back(x_t0->get(indices[coupling_id]));
        coupling_id++;
    }
    return local_xs_t0;
}

std::vector<std::vector<double>> getCurrentLocalSolutions(
    const CoupledSolutionsForStaggeredScheme& cpl_xs,
    const std::vector<std::vector<GlobalIndexType>>& indices)
{
    if (cpl_xs.coupled_xs.empty())
    {
        return {};
    }

    const auto number_of_coupled_solutions = cpl_xs.coupled_xs.size();
    std::vector<std::vector<double>> local_xs_t1;
    local_xs_t1.reserve(number_of_coupled_solutions);

    int coupling_id = 0;
    for (auto const* x_t1 : cpl_xs.coupled_xs)
    {
        local_xs_t1.emplace_back(x_t1->get(indices[coupling_id]));
        coupling_id++;
    }
    return local_xs_t1;
}

}  // namespace ProcessLib
