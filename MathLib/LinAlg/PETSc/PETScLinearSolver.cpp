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
    KSPSetOptionsPrefix(_solver, prefix.c_str());    
    KSPSetFromOptions(_solver);  // set running time option
}

void PETScLinearSolver::solve(const PETScVector &b, PETScVector &x)
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

    if(reason == KSP_DIVERGED_INDEFINITE_PC)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\nDivergence because of indefinite preconditioner;\n");
        PetscPrintf(PETSC_COMM_WORLD,"Run the executable again but with -pc_factor_shift_positive_definite option.\n");
    }
    else if(reason < 0)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\nOther kind of divergence: this should not happen.\n");
    }
    else
    {
        const char *slv_type;
        const char *prc_type;
        KSPGetType(_solver, &slv_type);
        PCGetType(_pc, &prc_type);

        PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
        PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner",
                    slv_type, prc_type);

        int its;
        KSPGetIterationNumber(_solver, &its);
        PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n", its);
        PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");

#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

} //end of namespace

