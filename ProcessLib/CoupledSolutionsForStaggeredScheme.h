/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CoupledSolutionsForStaggeredScheme.h
 *
 * Created on November 7, 2016, 12:14 PM
 */

#pragma once

#include <functional>
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
    CoupledSolutionsForStaggeredScheme(
        std::vector<std::reference_wrapper<GlobalVector const>> const&
            coupled_xs_,
        const double dt_, const int process_id_);

    /// References to the current solutions of the coupled processes.
    std::vector<std::reference_wrapper<GlobalVector const>> const& coupled_xs;

    /// Pointers to the vector of the solutions of the previous time step.
    std::vector<GlobalVector*> coupled_xs_t0;

    const double dt;  ///< Time step size.
    const int process_id;
};

/**
 *  A struct to keep the references to the local element solutions  of the
 *  current and previous time step solutions of the equations of the coupled
 *  processes.
 *
 *  During the global assembly loop, an instance of this struct is created for
 *  each element and it is then passed to local assemblers.
 */
struct LocalCoupledSolutions
{
    LocalCoupledSolutions(const double dt_, const int process_id_,
                          std::vector<std::vector<double>>&& local_coupled_xs0_,
                          std::vector<std::vector<double>>&& local_coupled_xs_)
        : dt(dt_),
          process_id(process_id_),
          local_coupled_xs0(std::move(local_coupled_xs0_)),
          local_coupled_xs(std::move(local_coupled_xs_))
    {
    }

    const double dt;  ///< Time step size.
    const int process_id;

    /// Local solutions of the previous time step.
    std::vector<std::vector<double>> const local_coupled_xs0;
    /// Local solutions of the current time step.
    std::vector<std::vector<double>> const local_coupled_xs;
};

std::vector<std::vector<double>> getPreviousLocalSolutions(
    const CoupledSolutionsForStaggeredScheme& cpl_xs,
    const std::vector<GlobalIndexType>& indices);

std::vector<std::vector<double>> getCurrentLocalSolutions(
    const CoupledSolutionsForStaggeredScheme& cpl_xs,
    const std::vector<GlobalIndexType>& indices);
}  // end of ProcessLib
