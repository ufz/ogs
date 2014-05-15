/*!
   \file  PETScLinearSolver.cpp
   \brief Definition of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Mar 2014

   \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScLinearSolver.h"

namespace MathLib
{
PETScLinearSolver::PETScLinearSolver(PETScMatrix &A, const std::string &prefix) : _A(A)
{
    KSPCreate(PETSC_COMM_WORLD, &_solver);

    KSPGetPC(_solver, &_pc);

    //
    if ( !prefix.empty() )
    {
        KSPSetOptionsPrefix(_solver, prefix.c_str());
    }

    KSPSetFromOptions(_solver);  // set running time option
}

bool PETScLinearSolver::solve(const PETScVector &b, PETScVector &x)
{
// define TEST_MEM_PETSC
#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    KSPSetOperators(_solver, _A.getRawMatrix(), _A.getRawMatrix(), DIFFERENT_NONZERO_PATTERN);

    KSPSolve(_solver, b.getData(), x.getData());

    KSPConvergedReason reason;
    KSPGetConvergedReason(_solver, &reason);

    bool converged = true;
    if(reason > 0)
    {
        const char *ksp_type;
        const char *pc_type;
        KSPGetType(_solver, &ksp_type);
        PCGetType(_pc, &pc_type);

        PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
        PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner",
                    ksp_type, pc_type);

        PetscInt its;
        KSPGetIterationNumber(_solver, &its);
        PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n", its);
        PetscPrintf(PETSC_COMM_WORLD,"\n================================================\n");
    }
    else if(reason == KSP_DIVERGED_ITS)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\nWarning: maximum number of iterations reached.\n");
    }
    else
    {
        converged = false;
        if(reason == KSP_DIVERGED_INDEFINITE_PC)
        {
            PetscPrintf(PETSC_COMM_WORLD, "\nDivergence because of indefinite preconditioner,");
            PetscPrintf(PETSC_COMM_WORLD, "\nTry to run again with -pc_factor_shift_positive_definite option.\n");
        }
        else if(reason == KSP_DIVERGED_BREAKDOWN_BICG)
        {
            PetscPrintf(PETSC_COMM_WORLD, "\nKSPBICG method was detected so the method could not continue to enlarge the Krylov space.");
            PetscPrintf(PETSC_COMM_WORLD, "\nTry to run again with another solver.\n");
        }
        else if(reason == KSP_DIVERGED_NONSYMMETRIC)
        {
            PetscPrintf(PETSC_COMM_WORLD, "\nMatrix or preconditioner is unsymmetric but KSP requires symmetric.\n");
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD, "\nDivergence detected, use command option -ksp_monitor or -log_summary to check the details.\n");
        }
    }

#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif

    return converged;
}

} //end of namespace

