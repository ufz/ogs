/*!
   \file
   \brief Definition of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Mar 2014

   \copyright
   Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScLinearSolver.h"

#include "BaseLib/RunTime.h"

namespace MathLib
{
PETScLinearSolver::PETScLinearSolver(std::string const& prefix,
                                     std::string const& petsc_options)
{
#if PETSC_VERSION_LT(3, 7, 0)
    PetscOptionsInsertString(petsc_options.c_str());
#else
    PetscOptionsInsertString(nullptr, petsc_options.c_str());
#endif

    KSPCreate(PETSC_COMM_WORLD, &solver_);

    KSPGetPC(solver_, &pc_);

    if (!prefix.empty())
    {
        KSPSetOptionsPrefix(solver_, prefix.c_str());
    }

    KSPSetInitialGuessNonzero(solver_, PETSC_TRUE);
    KSPSetFromOptions(solver_);  // set run-time options
}

bool PETScLinearSolver::solve(PETScMatrix& A, PETScVector& b, PETScVector& x)
{
    BaseLib::RunTime wtimer;
    wtimer.start();

// define TEST_MEM_PETSC
#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    KSPNormType norm_type;
    KSPGetNormType(solver_, &norm_type);
    const char* ksp_type;
    const char* pc_type;
    KSPGetType(solver_, &ksp_type);
    PCGetType(pc_, &pc_type);

    PetscPrintf(PETSC_COMM_WORLD,
                "\n================================================");
    PetscPrintf(PETSC_COMM_WORLD,
                "\nLinear solver %s with %s preconditioner using %s", ksp_type,
                pc_type, KSPNormTypes[norm_type]);

    KSPSetOperators(solver_, A.getRawMatrix(), A.getRawMatrix());

    KSPSolve(solver_, b.getRawVector(), x.getRawVector());

    KSPConvergedReason reason;
    KSPGetConvergedReason(solver_, &reason);

    bool converged = true;
    if (reason > 0)
    {
        PetscInt its;
        KSPGetIterationNumber(solver_, &its);
        PetscPrintf(PETSC_COMM_WORLD, "\nconverged in %d iterations", its);
        switch (reason)
        {
            case KSP_CONVERGED_RTOL:
                PetscPrintf(PETSC_COMM_WORLD,
                            " (relative convergence criterion fulfilled).");
                break;
            case KSP_CONVERGED_ATOL:
                PetscPrintf(PETSC_COMM_WORLD,
                            " (absolute convergence criterion fulfilled).");
                break;
            default:
                PetscPrintf(PETSC_COMM_WORLD, ".");
        }

        PetscPrintf(PETSC_COMM_WORLD,
                    "\n================================================\n");
    }
    else if (reason == KSP_DIVERGED_ITS)
    {
        PetscPrintf(PETSC_COMM_WORLD,
                    "\nWarning: maximum number of iterations reached.\n");
    }
    else
    {
        converged = false;
        if (reason == KSP_DIVERGED_INDEFINITE_PC)
        {
            PetscPrintf(PETSC_COMM_WORLD,
                        "\nDivergence because of indefinite preconditioner,");
            PetscPrintf(PETSC_COMM_WORLD,
                        "\nTry to run again with "
                        "-pc_factor_shift_positive_definite option.\n");
        }
        else if (reason == KSP_DIVERGED_BREAKDOWN_BICG)
        {
            PetscPrintf(PETSC_COMM_WORLD,
                        "\nKSPBICG method was detected so the method could not "
                        "continue to enlarge the Krylov space.");
            PetscPrintf(PETSC_COMM_WORLD,
                        "\nTry to run again with another solver.\n");
        }
        else if (reason == KSP_DIVERGED_NONSYMMETRIC)
        {
            PetscPrintf(PETSC_COMM_WORLD,
                        "\nMatrix or preconditioner is unsymmetric but KSP "
                        "requires symmetric.\n");
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD,
                        "\nDivergence detected, use command option "
                        "-ksp_monitor or -log_summary to check the details.\n");
        }
    }

#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(
        PETSC_COMM_WORLD,
        "###Memory usage by solver. Before: %f After: %f Increase: %d\n", mem1,
        mem2, (int)(mem2 - mem1));
#endif

    elapsed_ctime_ += wtimer.elapsed();

    return converged;
}

}  // namespace MathLib
