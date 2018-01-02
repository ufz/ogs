/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "LisLinearSolver.h"

#include <logog/include/logog.hpp>

#include "LisCheck.h"
#include "LisMatrix.h"
#include "LisVector.h"

namespace MathLib
{

LisLinearSolver::LisLinearSolver(
                    const std::string /*solver_name*/,
                    const BaseLib::ConfigTree* const option)
: _lis_option(option)
{
}

bool LisLinearSolver::solve(LisMatrix &A, LisVector &b, LisVector &x)
{
    finalizeMatrixAssembly(A);

    INFO("------------------------------------------------------------------");
    INFO("*** LIS solver computation");

    // Create solver
    LIS_SOLVER solver;
    int ierr = lis_solver_create(&solver);
    if (!checkLisError(ierr))
        return false;

    lis_solver_set_option(
        const_cast<char*>(_lis_option._option_string.c_str()), solver);
#ifdef _OPENMP
    INFO("-> number of threads: %i", (int) omp_get_max_threads());
#endif
    {
        int precon;
        ierr = lis_solver_get_precon(solver, &precon);
        INFO("-> precon: %i", precon);
    }
    {
        int slv;
        ierr = lis_solver_get_solver(solver, &slv);
        INFO("-> solver: %i", slv);
    }

    // solve
    INFO("-> solve");
    ierr = lis_solve(A.getRawMatrix(), b.getRawVector(), x.getRawVector(), solver);
    if (!checkLisError(ierr))
        return false;

    LIS_INT linear_solver_status;
    ierr = lis_solver_get_status(solver, &linear_solver_status);
    if (!checkLisError(ierr))
        return false;

    INFO("-> status: %d", linear_solver_status);

    {
        int iter = 0;
        ierr = lis_solver_get_iter(solver, &iter);
        if (!checkLisError(ierr))
            return false;

        INFO("-> iteration: %d", iter);
    }
    {
        double resid = 0.0;
        ierr = lis_solver_get_residualnorm(solver, &resid);
        if (!checkLisError(ierr))
            return false;
        INFO("-> residual: %g", resid);
    }
    {
        double time, itime, ptime, p_ctime, p_itime;
        ierr = lis_solver_get_timeex(solver, &time, &itime,
                                     &ptime, &p_ctime, &p_itime);
        if (!checkLisError(ierr))
            return false;
        INFO("-> time total           (s): %g", time);
        INFO("-> time iterations      (s): %g", itime);
        INFO("-> time preconditioning (s): %g", ptime);
        INFO("-> time precond. create (s): %g", p_ctime);
        INFO("-> time precond. iter   (s): %g", p_itime);
    }

    // Clear solver
    ierr = lis_solver_destroy(solver);
    if (!checkLisError(ierr))
        return false;
    INFO("------------------------------------------------------------------");

    return linear_solver_status == LIS_SUCCESS;
}

} //MathLib
