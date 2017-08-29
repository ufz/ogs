/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   StaggeredCouplingTerm.h
 *
 * Created on November 7, 2016, 12:14 PM
 */

#pragma once

#include <unordered_map>
#include <typeindex>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
class Process;

/**
 *  A struct to keep the references of the coupled processes and the references
 *  of the current solutions of the equations of the coupled processes.
 *
 *  During staggered coupling iteration, an instance of this struct is created
 *  and passed through interfaces to global and local assemblers for each
 *  process.
 */
struct StaggeredCouplingTerm
{
    StaggeredCouplingTerm(
        std::unordered_map<std::type_index, Process const&> const&
            coupled_processes_,
        std::unordered_map<std::type_index, GlobalVector const&> const&
            coupled_xs_,
        const double dt_);

    /// References to the coupled processes are distinguished by the keys of
    /// process types.
    std::unordered_map<std::type_index, Process const&> const&
        coupled_processes;

    /// References to the current solutions of the coupled processes.
    /// The coupled solutions are distinguished by the keys of process types.
    std::unordered_map<std::type_index, GlobalVector const&> const& coupled_xs;

    const double dt;  ///< Time step size.
};

/**
 *  A struct to keep the references of the coupled processes, and local element
 *  solutions  of the current and previous time step solutions of the equations
 *  of the coupled processes.
 *
 *  During the global assembly loop, an instance of this struct is created for
 *  each element and it is then passed to local assemblers.
 */
struct LocalCouplingTerm
{
    LocalCouplingTerm(
        const double dt_,
        std::unordered_map<std::type_index, Process const&> const&
            coupled_processes_,
        std::unordered_map<std::type_index, const std::vector<double>>&&
            local_coupled_xs0_,
        std::unordered_map<std::type_index, const std::vector<double>>&&
            local_coupled_xs_)
        : dt(dt_),
          coupled_processes(coupled_processes_),
          local_coupled_xs0(std::move(local_coupled_xs0_)),
          local_coupled_xs(std::move(local_coupled_xs_))
    {
    }

    const double dt;  ///< Time step size.

    /// References to the coupled processes are distinguished by the keys of
    /// process types.
    std::unordered_map<std::type_index, Process const&> const&
        coupled_processes;

    /// Local solutions of the previous time step.
    std::unordered_map<std::type_index, const std::vector<double>> const
        local_coupled_xs0;
    /// Local solutions of the current time step.
    std::unordered_map<std::type_index, const std::vector<double>> const
        local_coupled_xs;
};

std::unordered_map<std::type_index, const std::vector<double>>
getCurrentLocalSolutionsOfCoupledProcesses(
    const std::unordered_map<std::type_index, GlobalVector const&>&
        global_coupled_xs,
    const std::vector<GlobalIndexType>& indices);

}  // end of ProcessLib
