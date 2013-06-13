/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisLinearSolver.h"

#include <string>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "logog/include/logog.hpp"

#include "LisCheck.h"

namespace MathLib
{

void LisLinearSolver::setOption(const boost::property_tree::ptree &option)
{
    using boost::property_tree::ptree;

    boost::optional<ptree> ptSolver = option.get_child("LinearSolver");
    if (!ptSolver)
        return;

    boost::optional<std::string> solver_type = ptSolver->get_optional<std::string>("solver_type");
    if (solver_type) {
        _option.solver_type = _option.getSolverType(*solver_type);
    }
    boost::optional<std::string> precon_type = ptSolver->get_optional<std::string>("precon_type");
    if (precon_type) {
        _option.precon_type = _option.getPreconType(*precon_type);
    }
    boost::optional<std::string> matrix_type = ptSolver->get_optional<std::string>("matrix_type");
    if (matrix_type) {
        _option.matrix_type = _option.getMatrixType(*matrix_type);
    }
    boost::optional<double> error_tolerance = ptSolver->get_optional<double>("error_tolerance");
    if (error_tolerance) {
        _option.error_tolerance = *error_tolerance;
    }
    boost::optional<int> max_iteration_step = ptSolver->get_optional<int>("max_iteration_step");
    if (max_iteration_step) {
        _option.max_iterations = *max_iteration_step;
    }
}

void LisLinearSolver::solve(LisMatrix &A, LisVector &b, LisVector &x)
{
    if (!isMatrixAssembled(A)) {
        ERR("-> LisMatrix has not been assembled. LisLinearSolver::solve() is skipped.");
        return;
    }

    INFO("------------------------------------------------------------------");
    INFO("*** LIS solver computation");
#ifdef _OPENMP
    INFO("-> max number of threads = %d", omp_get_num_procs());
    INFO("-> number of threads = %d", omp_get_max_threads());
#endif

    int ierr = 0;

    // configure option
    std::string solver_options;
    if (_option.solver_precon_arg.empty()) {
        std::stringstream ss;
        ss << "-i " << static_cast<int>(_option.solver_type);
        ss << " -p " << static_cast<int>(_option.precon_type);
        if (!_option.extra_arg.empty())
            ss << " " << _option.extra_arg;
        solver_options = ss.str();
    } else {
        solver_options = _option.solver_precon_arg;
    }
    std::string tol_option;
    {
        std::stringstream ss;
        ss << "-tol " << _option.error_tolerance;
        ss << " -maxiter " << _option.max_iterations;
        ss << " -initx_zeros 0"; //0: use given x as initial guess, 1: x0=0
#ifdef _OPENMP
        const int nthreads = omp_get_max_threads();
        ss << " -omp_num_threads " << nthreads;
#endif
        tol_option = ss.str();
    }

    // Create solver
    LIS_SOLVER solver;
    ierr = lis_solver_create(&solver); checkLisError(ierr);
    ierr = lis_solver_set_option(const_cast<char*>(solver_options.c_str()), solver); checkLisError(ierr);
    ierr = lis_solver_set_option(const_cast<char*>(tol_option.c_str()), solver); checkLisError(ierr);
    ierr = lis_solver_set_option(const_cast<char*>("-print mem"), solver); checkLisError(ierr);

    // solve
    INFO("-> solve");
    ierr = lis_solve(A.getRawMatrix(), b.getRawVector(), x.getRawVector(), solver); checkLisError(ierr);

    int iter = 0;
    double resid = 0.0;
    ierr = lis_solver_get_iters(solver, &iter); checkLisError(ierr);
    ierr = lis_solver_get_residualnorm(solver, &resid); checkLisError(ierr);
    INFO("\t iteration: %d/%ld\n", iter, _option.max_iterations);
    INFO("\t residuals: %e\n", resid);

    // Clear solver
    ierr = lis_solver_destroy(solver); checkLisError(ierr);
    INFO("------------------------------------------------------------------");
}

} //MathLib
