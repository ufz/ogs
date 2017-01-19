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

#include <map>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

#include "ProcessType.h"

namespace ProcessLib
{
class Process;
struct StaggeredCouplingTerm
{
    StaggeredCouplingTerm(
        std::map<ProcessType, Process const&> const& coupled_processes_,
        std::map<ProcessType, GlobalVector const*> const& coupled_xs_,
        const double dt_, const bool empty_ = false)
    : coupled_processes(coupled_processes_), coupled_xs(coupled_xs_),
      dt(dt_), empty(empty_)
    {
    }

    std::map<ProcessType, Process const&> const& coupled_processes;
    std::map<ProcessType, GlobalVector const*> const& coupled_xs;
    const double dt; ///< Time step size.
    const bool empty;
};

struct LocalCouplingTerm
{
    LocalCouplingTerm(const double dt_,
        std::map<ProcessType, Process const&> const& coupled_processes_,
        std::map<ProcessType, const std::vector<double>>&& local_coupled_xs0_,
        std::map<ProcessType, const std::vector<double>>&& local_coupled_xs_)
    : dt(dt_), coupled_processes(coupled_processes_),
      local_coupled_xs0(std::move(local_coupled_xs0_)),
      local_coupled_xs(std::move(local_coupled_xs_))
    {
    }

    const double dt; ///< Time step size.

    std::map<ProcessType, Process const&> const& coupled_processes;
    /// Local solutions of the previous time step.
    std::map<ProcessType, const std::vector<double>> const local_coupled_xs0;
    /// Local solutions of the current time step.
    std::map<ProcessType, const std::vector<double>> const local_coupled_xs;
};

const StaggeredCouplingTerm createVoidStaggeredCouplingTerm();

} // end of ProcessLib

