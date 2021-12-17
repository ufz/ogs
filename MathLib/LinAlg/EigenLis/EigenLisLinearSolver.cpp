/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLisLinearSolver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#include "MathLib/LinAlg/Lis/LisVector.h"

namespace MathLib
{
EigenLisLinearSolver::EigenLisLinearSolver(
    const std::string /*solver_name*/, BaseLib::ConfigTree const* const option)
    : lis_option_(option)
{
}

EigenLisLinearSolver::EigenLisLinearSolver(std::string const& /*solver_name*/,
                                           std::string const& lis_options)
    : lis_option_(nullptr), lis_options_(lis_options)
{
}

bool EigenLisLinearSolver::solve(EigenMatrix& A_, EigenVector& b_,
                                 EigenVector& x_)
{
    static_assert(EigenMatrix::RawMatrixType::IsRowMajor,
                  "Sparse matrix is required to be in row major storage.");
    auto& A = A_.getRawMatrix();
    auto& b = b_.getRawVector();
    auto& x = x_.getRawVector();

    if (!A.isCompressed())
    {
        A.makeCompressed();
    }
    int nnz = A.nonZeros();
    int* ptr = A.outerIndexPtr();
    int* col = A.innerIndexPtr();
    double* data = A.valuePtr();
    LisMatrix lisA(A_.getNumberOfRows(), nnz, ptr, col, data);
    LisVector lisb(b.rows(), b.data());
    LisVector lisx(x.rows(), x.data());

    bool const status = solve(lisA, lisb, lisx);

    for (std::size_t i = 0; i < lisx.size(); i++)
    {
        x[i] = lisx[i];
    }

    return status;
}

bool EigenLisLinearSolver::solve(LisMatrix& A, LisVector& b, LisVector& x)
{
    finalizeMatrixAssembly(A);

    INFO("------------------------------------------------------------------");
    INFO("*** LIS solver computation");

    // Create solver
    LIS_SOLVER solver;
    int ierr = lis_solver_create(&solver);
    if (!checkLisError(ierr))
    {
        return false;
    }
    lis_solver_set_option(const_cast<char*>(lis_options_.c_str()), solver);
#ifdef _OPENMP
    INFO("-> number of threads: {:d}", (int)omp_get_max_threads());
#endif
    {
        int precon;
        ierr = lis_solver_get_precon(solver, &precon);
        if (!checkLisError(ierr))
        {
            return false;
        }
        INFO("-> precon: {:d}", precon);
    }
    {
        int slv;
        ierr = lis_solver_get_solver(solver, &slv);
        if (!checkLisError(ierr))
        {
            return false;
        }
        INFO("-> solver: {:d}", slv);
    }

    // solve
    INFO("-> solve");
    ierr =
        lis_solve(A.getRawMatrix(), b.getRawVector(), x.getRawVector(), solver);
    if (!checkLisError(ierr))
    {
        return false;
    }

    LIS_INT linear_solver_status;
    ierr = lis_solver_get_status(solver, &linear_solver_status);
    if (!checkLisError(ierr))
    {
        return false;
    }

    INFO("-> status: {:d}", linear_solver_status);

    {
        int iter = 0;
        ierr = lis_solver_get_iter(solver, &iter);
        if (!checkLisError(ierr))
        {
            return false;
        }

        INFO("-> iteration: {:d}", iter);
    }
    {
        double resid = 0.0;
        ierr = lis_solver_get_residualnorm(solver, &resid);
        if (!checkLisError(ierr))
        {
            return false;
        }
        INFO("-> residual: {:g}", resid);
    }
    {
        double time, itime, ptime, p_ctime, p_itime;
        ierr = lis_solver_get_timeex(solver, &time, &itime, &ptime, &p_ctime,
                                     &p_itime);
        if (!checkLisError(ierr))
        {
            return false;
        }
        INFO("-> time total           (s): {:g}", time);
        INFO("-> time iterations      (s): {:g}", itime);
        INFO("-> time preconditioning (s): {:g}", ptime);
        INFO("-> time precond. create (s): {:g}", p_ctime);
        INFO("-> time precond. iter   (s): {:g}", p_itime);
    }

    // Clear solver
    ierr = lis_solver_destroy(solver);
    if (!checkLisError(ierr))
    {
        return false;
    }
    INFO("------------------------------------------------------------------");

    return linear_solver_status == LIS_SUCCESS;
}

}  // namespace MathLib
