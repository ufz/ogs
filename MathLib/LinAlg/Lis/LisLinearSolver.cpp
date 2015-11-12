/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisLinearSolver.h"

#include "logog/include/logog.hpp"

#include "LisCheck.h"

namespace MathLib
{

LisLinearSolver::LisLinearSolver(LisMatrix &A,
                    const std::string /*solver_name*/,
                    BaseLib::ConfigTree const*const option)
: _A(A)
{
    if (option)
    {
        _option.addOptions(*option);

#ifndef NDEBUG
        _option.printInfo();
#endif
    }
}

void LisLinearSolver::solve(LisVector &b, LisVector &x)
{
    finalizeMatrixAssembly(_A);

    INFO("------------------------------------------------------------------");
    INFO("*** LIS solver computation");

    // Create solver
    LIS_SOLVER solver;
    int ierr = lis_solver_create(&solver);
    checkLisError(ierr);

    {
        std::string opt;
        for (auto const& it : _option.settings)
        {
            opt = it.first + " " + it.second;

            ierr = lis_solver_set_option(const_cast<char*>(opt.c_str()), solver);
            checkLisError(ierr);
        }
    }
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
    ierr = lis_solve(_A.getRawMatrix(), b.getRawVector(), x.getRawVector(), solver);
    checkLisError(ierr);

    {
        int iter = 0;
        ierr = lis_solver_get_iter(solver, &iter);
        checkLisError(ierr);

        std::string max_iter = _option.settings["-maxiter"];
        if (max_iter.empty()) max_iter = "--";
        INFO("-> iteration: %d/%s", iter, max_iter.c_str());
    }
    {
        double resid = 0.0;
        ierr = lis_solver_get_residualnorm(solver, &resid);
        checkLisError(ierr);
        INFO("-> residual: %g", resid);
    }
    {
        double time, itime, ptime, p_ctime, p_itime;
        ierr = lis_solver_get_timeex(solver, &time, &itime,
                                     &ptime, &p_ctime, &p_itime);
        checkLisError(ierr);
        INFO("-> time total           (s): %g", time);
        INFO("-> time iterations      (s): %g", itime);
        INFO("-> time preconditioning (s): %g", ptime);
        INFO("-> time precond. create (s): %g", p_ctime);
        INFO("-> time precond. iter   (s): %g", p_itime);
    }

    // Clear solver
    ierr = lis_solver_destroy(solver);
    checkLisError(ierr);
    INFO("------------------------------------------------------------------");
}

} //MathLib
