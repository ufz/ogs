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

#include <numeric>

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

std::vector<double> getPreviousLocalSolutions(
    const CoupledSolutionsForStaggeredScheme& cpl_xs,
    const std::vector<std::vector<GlobalIndexType>>& indices)
{
    if (cpl_xs.coupled_xs_t0.empty())
    {
        return {};
    }

    std::size_t const local_solutions_size = std::accumulate(
        cbegin(indices),
        cend(indices),
        std::size_t(0),
        [](GlobalIndexType const size,
           std::vector<GlobalIndexType> const& process_indices) {
            return size + process_indices.size();
        });
    std::vector<double> local_xs_t0;
    local_xs_t0.reserve(local_solutions_size);

    int process_id = 0;
    for (auto const& x_t0 : cpl_xs.coupled_xs_t0)
    {
        auto const& values = x_t0->get(indices[process_id]);
        local_xs_t0.insert(cend(local_xs_t0), cbegin(values), cend(values));
        process_id++;
    }
    return local_xs_t0;
}

std::vector<double> getCurrentLocalSolutions(
    const CoupledSolutionsForStaggeredScheme& cpl_xs,
    const std::vector<std::vector<GlobalIndexType>>& indices)
{
    if (cpl_xs.coupled_xs.empty())
    {
        return {};
    }

    std::size_t const local_solutions_size = std::accumulate(
        cbegin(indices),
        cend(indices),
        std::size_t(0),
        [](GlobalIndexType const size,
           std::vector<GlobalIndexType> const& process_indices) {
            return size + process_indices.size();
        });
    std::vector<double> local_xs_t1;
    local_xs_t1.reserve(local_solutions_size);

    int process_id = 0;
    for (auto const& x_t1 : cpl_xs.coupled_xs)
    {
        auto const& values = x_t1->get(indices[process_id]);
        local_xs_t1.insert(cend(local_xs_t1), cbegin(values), cend(values));
        process_id++;
    }
    return local_xs_t1;
}

}  // namespace ProcessLib
