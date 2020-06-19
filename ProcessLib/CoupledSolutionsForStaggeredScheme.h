/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 7, 2016, 12:14 PM
 */

#pragma once

#include <utility>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
/**
 *  A struct to keep the references of the current solutions of the equations of
 *  the coupled processes.
 *
 *  During staggered coupling iteration, an instance of this struct is created
 *  and passed through interfaces to global and local assemblers for each
 *  process.
 */
struct CoupledSolutionsForStaggeredScheme
{
    explicit CoupledSolutionsForStaggeredScheme(
        std::vector<GlobalVector*> const& coupled_xs_);

    /// References to the current solutions of the coupled processes.
    std::vector<GlobalVector*> const& coupled_xs;

    /// Pointers to the vector of the solutions of the previous time step.
    std::vector<GlobalVector*> coupled_xs_t0;
};

/**
    LocalCoupledSolutions(std::vector<double>&& local_coupled_xs0_)
 * Fetch the nodal solutions of all coupled processes from the given vector of
 * global solutions for each process into a flat vector.
 */
std::vector<double> getCoupledLocalSolutions(
    std::vector<GlobalVector*> const& global_solutions,
    std::vector<std::vector<GlobalIndexType>> const& indices);

}  // namespace ProcessLib
