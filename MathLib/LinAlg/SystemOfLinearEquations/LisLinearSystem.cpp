/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file LisLinearSystem.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "LisLinearSystem.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include "logog/include/logog.hpp"

namespace MathLib
{

bool LisLinearSystem::checkError(int err)
{
    bool ok = (err == LIS_SUCCESS);
    if (!ok) {
        ERR("***ERROR: Lis error code = %d", err);
    }
    return ok;
}

LisLinearSystem::LisLinearSystem(std::size_t dimension, RowMajorSparsity* /*sp*/)
: _dim(dimension), _is_initialized(false), _max_diag_coeff(0)
{
    setZero();
}

LisLinearSystem::~LisLinearSystem()
{
    lis_matrix_destroy(_AA);
    lis_vector_destroy(_bb);
    lis_vector_destroy(_xx);
}

void LisLinearSystem::setOption(const boost::property_tree::ptree &option)
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

void LisLinearSystem::setZero()
{
    int ierr = 0;
    // create a matrix
    if (_is_initialized) {
        // Matrix has to be removed and created every time because Lis doesn't provide a
        // function to set matrix entries to zero
        ierr = lis_matrix_destroy(_AA);
    }
#ifndef USE_MPI
    ierr = lis_matrix_create(0, &_AA);
    ierr = lis_matrix_set_size(_AA, 0, getDimension());
#else
    lis_matrix_create(MPI_COMM_WORLD, &_AA);
    lis_matrix_set_size(_AA, dimension, 0);
//    lis_matrix_get_range(A,&is,&ie);
#endif
    checkError(ierr);

    // crate or zero RHS, x
    if (!_is_initialized) {
        ierr = lis_vector_duplicate(_AA, &_bb);
        ierr = lis_vector_duplicate(_AA, &_xx);
        _is_initialized = true;
    } else {
        ierr = lis_vector_set_all(0.0, _bb);
        ierr = lis_vector_set_all(0.0, _xx);
    }
    checkError(ierr);

    _max_diag_coeff = .0;
}

void LisLinearSystem::applyKnownSolution()
{
    //Use penalty parameter
    const double penalty_scaling = 1e+10;
    const double penalty = _max_diag_coeff * penalty_scaling;
    INFO("-> max. absolute value of diagonal entries = %e", _max_diag_coeff);
    INFO("-> penalty scaling = %e", penalty_scaling);
    const std::size_t n_bc = _vec_knownX_id.size();
    for (std::size_t i_bc=0; i_bc<n_bc; i_bc++) {
        const std::size_t rowId = _vec_knownX_id[i_bc];
        const double x = _vec_knownX_x[i_bc];

        //A(k, k) = epsilon
        setMatEntry(rowId, rowId, penalty);
        //b(k) = x*epsilon
        setRHSVec(rowId, x*penalty);
    }
}

void LisLinearSystem::solve()
{
    INFO("------------------------------------------------------------------");
    INFO("*** LIS solver computation");
#ifdef _OPENMP
    INFO("-> max number of threads = %d", omp_get_num_procs());
    INFO("-> number of threads = %d", omp_get_max_threads());
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif

    int ierr = 0;

    applyKnownSolution();
    // assemble a matrix
    ierr = lis_matrix_set_type(_AA, static_cast<int>(_option.matrix_type));
    ierr = lis_matrix_assemble(_AA);
    checkError(ierr);

    // configure option
    const std::size_t MAX_ZEILE = 512;
    char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];
    if (_option.solver_precon_arg.length()==0) {
        sprintf(solver_options, "-i %d -p %d %s ", static_cast<int>(_option.solver_type), static_cast<int>(_option.precon_type), _option.extra_arg.c_str());
    } else {
        sprintf(solver_options, "%s ", _option.solver_precon_arg.c_str());
    }
    sprintf(tol_option, "-tol %e -maxiter %ld -omp_num_threads %d -initx_zeros 0", _option.error_tolerance, _option.max_iterations, nthreads);

    // Create solver
    LIS_SOLVER solver;
    ierr = lis_solver_create(&solver);
    ierr = lis_solver_set_option(solver_options, solver);
    ierr = lis_solver_set_option(tol_option, solver);
    ierr = lis_solver_set_option((char*)"-print mem", solver);
    
    // solve
    INFO("-> solve");
    ierr = lis_solve(_AA, _bb, _xx, solver);
    checkError(ierr);
    //lis_output(_AA, _bb, _xx, LIS_FMT_MM, "/home/localadmin/tasks/20120814_ogs6test/matrix1.txt");

    int iter = 0;
    double resid = 0.0;
    ierr = lis_solver_get_iters(solver, &iter);
    ierr = lis_solver_get_residualnorm(solver, &resid);
    printf("\t iteration: %d/%ld\n", iter, _option.max_iterations);
    printf("\t residuals: %e\n", resid);

    // Clear solver
    lis_solver_destroy(solver);
    std::cout << "------------------------------------------------------------------" << std::endl;
}

void LisLinearSystem::printout(std::ostream &os) const
{
    os << "#A=" << std::endl;
    os << "#x=" << std::endl;
    lis_vector_print(_xx);
    os << "#b=" << std::endl;
    lis_vector_print(_bb);
}

} //MathLib
