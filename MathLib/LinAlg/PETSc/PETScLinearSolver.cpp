/*!
   \file  PETScLinearSolve.cpp
   \brief Definition of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.


   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013


  \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/


#include "PETScLinearSolver.h"

#include<iostream>



namespace MathLib
{

using boost::property_tree::ptree;


PETScLinearSolver::PETScLinearSolver(const PetscInt size, PETScMatrix _matrix,  PETScVector _vec)
   :A(_matrix), b(_vec), x(_vec), lsolver(nullptr), prec(nullptr), global_x0(nullptr),  global_x1(nullptr)
{
   ltolerance = 1.e-10;
   time_elapsed = 0.0;
   m_size_loc = PETSC_DECIDE;
   m_size = size;


}


PETScLinearSolver::PETScLinearSolver(PETScMatrix &stiffness_matrix,
                                     boost::property_tree::ptree const*const option,
                                     PETScVector _vec)
   :A(stiffness_matrix), b(_vec), x(_vec), lsolver(nullptr), prec(nullptr), global_x0(nullptr),  global_x1(nullptr)
{
   ltolerance = 1.e-10;
   time_elapsed = 0.0;

   m_size = A.size();

   m_size_loc = PETSC_DECIDE;

   mpi_size = A.getMPI_Size();
   rank = A.getMPI_Rank();

   if (option)
      Config(*option);

   allocateMemory4TemoraryArrays(m_size);

}

PETScLinearSolver :: PETScLinearSolver (PETScMatrix &stiffness_matrix,
                                        PETScVector &rhs, PETScVector &unknowns)
   :A(stiffness_matrix), b(rhs), x(unknowns), lsolver(nullptr), prec(nullptr), global_x0(nullptr),  global_x1(nullptr)
{
   ltolerance = 1.e-10;
   time_elapsed = 0.0;

   m_size_loc = PETSC_DECIDE;
   m_size = A.size();

   mpi_size = A.getMPI_Size();
   rank = A.getMPI_Rank();

   x.set_rank_size(rank, mpi_size);
   b.set_rank_size(rank, mpi_size);


}



PETScLinearSolver:: ~PETScLinearSolver()
{
   if(lsolver) KSPDestroy(&lsolver);
   // if(prec) PCDestroy(&prec);

   if(global_x0 != nullptr)
      delete []  global_x0;
   if(global_x1 != nullptr)
      delete []  global_x1;


   global_x0 = nullptr;
   global_x1 = nullptr;

   PetscPrintf(PETSC_COMM_WORLD,"\n>>Number of Unknows: %d", m_size);
   PetscPrintf(PETSC_COMM_WORLD,"\n>>Elapsed time in linear solver: %fs", time_elapsed);
}


void PETScLinearSolver:: allocateMemory4TemoraryArrays(const PetscInt size)
{

   m_size = size;

   if(global_x0 == nullptr)
      global_x0 = new PetscScalar[m_size];

   if(global_x1 == nullptr)
      global_x1 = new PetscScalar[m_size];

}


void PETScLinearSolver:: releaseMemory4TemoraryArrays()
{

   if(global_x0 != nullptr)
      delete []  global_x0;
   if(global_x1 != nullptr)
      delete []  global_x1;


   global_x0 = nullptr;
   global_x1 = nullptr;

}



/*!
  \brief KSP and PC type

 KSPRICHARDSON "richardson"
 KSPCHEBYCHEV  "chebychev"
 KSPCG         "cg"
 KSPCGNE       "cgne"
 KSPNASH       "nash"
 KSPSTCG       "stcg"
 KSPGLTR       "gltr"
 KSPGMRES      "gmres"
 KSPFGMRES     "fgmres"
 KSPLGMRES     "lgmres"
 KSPDGMRES     "dgmres"
 KSPTCQMR      "tcqmr"
 KSPBCGS       "bcgs"
 KSPIBCGS        "ibcgs"
 KSPBCGSL        "bcgsl"
 KSPCGS        "cgs"
 KSPTFQMR      "tfqmr"
 KSPCR         "cr"
 KSPLSQR       "lsqr"
 KSPPREONLY    "preonly"
 KSPQCG        "qcg"
 KSPBICG       "bicg"
 KSPMINRES     "minres"
 KSPSYMMLQ     "symmlq"
 KSPLCD        "lcd"
 KSPPYTHON     "python"
 KSPBROYDEN    "broyden"
 KSPGCR        "gcr"
 KSPNGMRES     "ngmres"
 KSPSPECEST    "specest"

 PCNONE            "none"
 PCJACOBI          "jacobi"
 PCSOR             "sor"
 PCLU              "lu"
 PCSHELL           "shell"
 PCBJACOBI         "bjacobi"
 PCMG              "mg"
 PCEISENSTAT       "eisenstat"
 PCILU             "ilu"
 PCICC             "icc"
 PCASM             "asm"
 PCGASM            "gasm"
 PCKSP             "ksp"
 PCCOMPOSITE       "composite"
 PCREDUNDANT       "redundant"
 PCSPAI            "spai"
 PCNN              "nn"
 PCCHOLESKY        "cholesky"
 PCPBJACOBI        "pbjacobi"
 PCMAT             "mat"
 PCHYPRE           "hypre"
 PCPARMS           "parms"
 PCFIELDSPLIT      "fieldsplit"
 PCTFS             "tfs"
 PCML              "ml"
 PCPROMETHEUS      "prometheus"
 PCGALERKIN        "galerkin"
 PCEXOTIC          "exotic"
 PCHMPI            "hmpi"
 PCSUPPORTGRAPH    "supportgraph"
 PCASA             "asa"
 PCCP              "cp"
 PCBFBT            "bfbt"
 PCLSC             "lsc"
 PCPYTHON          "python"
 PCPFMG            "pfmg"
 PCSYSPFMG         "syspfmg"
 PCREDISTRIBUTE    "redistribute"
 PCSACUSP          "sacusp"
 PCSACUSPPOLY      "sacusppoly"
 PCBICGSTABCUSP    "bicgstabcusp"
 PCSVD             "svd"
 PCAINVCUSP        "ainvcusp"
 PCGAMG            "gamg"

*/
void PETScLinearSolver::Config(const boost::property_tree::ptree &option)

{

   int maxits = 0;


   boost::optional<ptree> ptSolver = option.get_child("LinearSolver");
   if (!ptSolver)
      return;

   boost::optional<std::string> solver_type = ptSolver->get_optional<std::string>("solver_type");
   if (solver_type)
   {
      sol_type  = *solver_type;
   }
   boost::optional<std::string> precon_type = ptSolver->get_optional<std::string>("precon_type");
   if (precon_type)
   {
      pc_type = *precon_type;
   }
   boost::optional<double> error_tolerance = ptSolver->get_optional<double>("error_tolerance");
   if (error_tolerance)
   {
      ltolerance = *error_tolerance;
   }
   boost::optional<int> max_iteration_step = ptSolver->get_optional<int>("max_iteration_step");
   if (max_iteration_step)
   {
      maxits = *max_iteration_step;
   }

   const KSPType lsol_type = sol_type.c_str();
   const PCType prec_type = pc_type.c_str();


   KSPCreate(PETSC_COMM_WORLD,&lsolver);
   KSPSetOperators(lsolver, A.getData(), A.getData(), DIFFERENT_NONZERO_PATTERN);
   KSPSetType(lsolver,lsol_type);

   KSPGetPC(lsolver, &prec);
   PCSetType(prec, prec_type); //  PCJACOBI); //PCNONE);
   KSPSetTolerances(lsolver,ltolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxits);
   KSPSetFromOptions(lsolver);



}


void PETScLinearSolver::Solver()
{

   //TEST
#ifdef TEST_MEM_PETSC
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif


   int its;
   PetscLogDouble v1,v2;
   KSPConvergedReason reason;

   PetscGetTime(&v1);

   KSPSolve(lsolver, b.getData(), x.getData());

   KSPGetConvergedReason(lsolver,&reason); //CHKERRQ(ierr);
   if (reason==KSP_DIVERGED_INDEFINITE_PC)
   {
      PetscPrintf(PETSC_COMM_WORLD,"\nDivergence because of indefinite preconditioner;\n");
      PetscPrintf(PETSC_COMM_WORLD,"Run the executable again but with -pc_factor_shift_positive_definite option.\n");
   }
   else if (reason<0)
   {
      PetscPrintf(PETSC_COMM_WORLD,"\nOther kind of divergence: this should not happen.\n");
   }
   else
   {
      const char *slv_type;
      const char *prc_type;
      KSPGetType(lsolver, &slv_type);
      PCGetType(prec, &prc_type);

      PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
      PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner",
                  slv_type, prc_type);
      KSPGetIterationNumber(lsolver,&its); //CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n",(int)its);
      PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
   }
   PetscPrintf(PETSC_COMM_WORLD,"\n");

   //VecAssemblyBegin(x);
   //VecAssemblyEnd(x);

   PetscGetTime(&v2);
   time_elapsed += v2-v1;


#ifdef TEST_MEM_PETSC
   //TEST
   PetscMemoryGetCurrentUsage(&mem2);
   PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}



void PETScLinearSolver::Solver(PETScVector &rhs, PETScVector &unknowns)
{
   b = rhs;
   x = unknowns;

   b.set_rank_size(rank, mpi_size);
   x.set_rank_size(rank, mpi_size);


   Solver();
}


void PETScLinearSolver::finalAssemble()
{
   A.finalAssemble(MAT_FINAL_ASSEMBLY);

   //if(full_fill_A_b_x)
   {
      b.finalAssemble();
   }

}

void PETScLinearSolver::addMatrixEntries(const PetscInt m, const PetscInt idxm[],                                         const PetscInt n, const PetscInt idxn[],
      const PetscScalar v[])
{

   A.addEntries(m, idxm, n, idxn, v);
}


void PETScLinearSolver::applyKnownSolutions(PetscInt ni,const PetscInt ix[], const PetscScalar y[])
{
   A.zeroRows_in_Matrix(ni, ix);
   A.finalAssemble();



   //if(full_fill_A_b_x)
   {
      x.setValues(ni, ix, y, INSERT_VALUES);
      b.setValues(ni, ix, y, INSERT_VALUES);

      x.finalAssemble();
      b.finalAssemble();
   }
}


void PETScLinearSolver::mappingSolution()
{
   x.getGlobalEntries(global_x0, global_x1);
}



PetscScalar *PETScLinearSolver::getGlobalSolution() const
{

   return global_x1;

}

PetscReal PETScLinearSolver::getNormRHS(NormType  nmtype)
{

   return b.getNorm(nmtype);
}
/*!
    Get norm of x
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i|
    06.2012. WW
*/
PetscReal PETScLinearSolver::getNormUnknowns(NormType  nmtype)
{
   return x.getNorm(nmtype);
}


void PETScLinearSolver::initializeMatVec( )
{
   A.setZero();

   //if(full_fill_A_b_x)
   {
      b.setZero();
      x.setZero();
   }
}



void PETScLinearSolver::Viewer(std::string file_name)
{
   std::string fname = file_name + "_eqs_dump_";
   A.Viewer(fname + "matrix");

   // if(full_fill_A_b_x)
   {
      b.Viewer(fname + "rhs");
      x.Viewer(fname + "unknowns");
   }


#define  nEXIT_TEST
#ifdef EXIT_TEST
   if(lsolver) KSPDestroy(&lsolver);
   // if(prec) PCDestroy(&prec);
   if(global_x0)
      delete []  global_x0;
   if(global_x1)
      delete []  global_x1;

   global_x0 = nullptr;
   global_x1 = nullptr;

   PetscFinalize();
   exit(0);
#endif

}

} //end of namespace

