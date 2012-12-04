/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file LisLinearEquation.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "LisLinearEquation.h"

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace MathLib
{

void LisLinearEquation::initialize()
{
    int argc=0;
    char **argv;
    lis_initialize(&argc, &argv);
}

void LisLinearEquation::finalize()
{
    lis_finalize();
}

LisLinearEquation::LisLinearEquation(std::size_t length, RowMajorSparsity* sp)
: AbstractCRSLinearEquation<signed>(length, sp)
{
    initialize();
}

LisLinearEquation::~LisLinearEquation()
{
    lis_vector_destroy(bb);
    lis_vector_destroy(xx);
}

void LisLinearEquation::setOption(const BaseLib::Options &option)
{
    const BaseLib::Options *op = option.getSubGroup("LinearSolver");
    if (op==0) {
        return;
    }

    if (op->hasOption("solver_type"))
        _option.ls_method = _option.getSolverType(op->getOption("solver_type"));
    if (op->hasOption("precon_type"))
        _option.ls_precond = _option.getPreconType(op->getOption("precon_type"));
    if (op->hasOption("error_tolerance"))
        _option.ls_error_tolerance = op->getOptionAsNum<double>("error_tolerance");
    if (op->hasOption("max_iteration_step"))
        _option.ls_max_iterations = op->getOptionAsNum<int>("max_iteration_step");
}

void LisLinearEquation::solveEqs(CRSMatrix<double, signed> *A, double *b, double *x)
{
    long dimension = static_cast<long>(A->getNRows());

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "*** LIS solver computation" << std::endl;

//    std::cout << "A=" << std::endl;
//    A->printMat();
//    std::cout << "b=" << std::endl;
//    for (size_t i=0; i<A->getNRows(); i++)
//        std::cout << b[i] << " ";
//    std::cout << std::endl;


    int ierr = 0;
    // Creating a matrix.
    LIS_SOLVER solver;
    LIS_MATRIX AA;
#ifndef USE_MPI
    ierr = lis_matrix_create(0, &AA);
    ierr = lis_matrix_set_type(AA, LIS_MATRIX_CRS);
    ierr = lis_matrix_set_size(AA, 0, dimension);
#else
    lis_matrix_create(MPI_COMM_WORLD, &AA);
    lis_matrix_set_size(AA,dimension,0);
//    lis_matrix_get_size(AA, &_local_dim, &_global_dim);
#endif

    // Matrix solver and Precondition can be handled better way.
    const std::size_t MAX_ZEILE = 512;
    char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
//    omp_set_num_threads (nthreads);
#else
    int nthreads = 1;
#endif

    sprintf(solver_options, "-i %d -p %d %s", _option.ls_method, _option.ls_precond, _option.ls_extra_arg.c_str()); 
    sprintf(tol_option, "-tol %e -maxiter %ld -omp_num_threads %d -initx_zeros 0", _option.ls_error_tolerance, _option.ls_max_iterations, nthreads);

    ierr = lis_matrix_set_crs(A->getNNZ(), (int*)A->getRowPtrArray(), (int*)A->getColIdxArray(), (double*)A->getEntryArray(), AA);
    ierr = lis_matrix_assemble(AA);

    // Assemble the vector, b, x
    ierr = lis_vector_duplicate(AA, &bb);
    ierr = lis_vector_duplicate(AA, &xx);
    #pragma omp parallel for
    for (long i=0; i < dimension; ++i)
    {
        ierr = lis_vector_set_value(LIS_INS_VALUE, i, x[i], xx);
        ierr = lis_vector_set_value(LIS_INS_VALUE, i, b[i], bb);
    }

    // Create solver
    ierr = lis_solver_create(&solver);

    ierr = lis_solver_set_option(solver_options, solver);
    ierr = lis_solver_set_option(tol_option, solver);
    ierr = lis_solver_set_option("-print mem", solver);
    
    ierr = lis_solve(AA, bb, xx, solver);
    //lis_output(AA, bb, xx, LIS_FMT_MM, "/home/norihiro/work/task/20120814_ogs6test/deformation/matrix1.txt");

    int iter = 0;
    double resid = 0.0;
    ierr = lis_solver_get_iters(solver, &iter);
    ierr = lis_solver_get_residualnorm(solver, &resid);
    printf("\t iteration: %d/%ld\n", iter, _option.ls_max_iterations);
    printf("\t residuals: %e\n", resid);
    //    lis_vector_print(xx);
    //    lis_vector_print(bb);

    // Update the solution (answer) into the x vector
    #pragma omp parallel for
    for(long i=0; i<dimension; ++i)
    {
        lis_vector_get_value(xx,i,&(x[i]));
    }

    // Clear memory
    //lis_matrix_destroy(AA); // CRSMatrix will free this memory
    lis_solver_destroy(solver);
    std::cout << "------------------------------------------------------------------" << std::endl;
}

//void LisLinearEquation::gatherX(std::vector<double> &x)
//{
//    int local_n, global_n;
//    lis_vector_get_size(xx, &local_n, &global_n);
//    x.resize(global_n);
//    lis_vector_gather(xx, &x[0]);
//};

} //MathLib
