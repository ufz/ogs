/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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

std::vector<double> getCoupledLocalSolutions(
    std::vector<GlobalVector*> const& global_solutions,
    std::vector<std::vector<GlobalIndexType>> const& indices)
{
    if (global_solutions.empty())
    {
        return {};
    }

    std::size_t const local_solutions_size =
        std::accumulate(cbegin(indices),
                        cend(indices),
                        std::size_t(0),
                        [](GlobalIndexType const size,
                           std::vector<GlobalIndexType> const& process_indices)
                        { return size + process_indices.size(); });
    std::vector<double> local_solutions;
    local_solutions.reserve(local_solutions_size);

    int number_of_processes = static_cast<int>(global_solutions.size());
    for (int process_id = 0; process_id < number_of_processes; ++process_id)
    {
        auto values = global_solutions[process_id]->get(indices[process_id]);
        local_solutions.insert(cend(local_solutions),
                               std::make_move_iterator(begin(values)),
                               std::make_move_iterator(end(values)));
    }
    return local_solutions;
}
}  // namespace ProcessLib
