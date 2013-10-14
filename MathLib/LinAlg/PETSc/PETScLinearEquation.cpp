/*!
   \file  PETScLinearEquation.cpp
   \brief Definition of member functions of class PETScLinearEquation, which provides interfaces to
          matrix and solvers of PETSc routines.

   \author Wenqing Wang
   \version 
   \date Nov 2011 - Sep 2013


  \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/


#include "PETScLinearEquation.h"

#include<iostream>


namespace MathLib
{

using boost::property_tree::ptree;

 PETScLinearEquation :: PETScLinearEquation ()
   :lsolver(NULL), prec(NULL), global_x0(NULL),  global_x1(NULL), global_buff(NULL) 
 {
   ltolerance = 1.e-10;
   time_elapsed = 0.0;   
   d_nz = 10; 
   o_nz = 10; 	
   nz = 10;   
   m_size_loc = PETSC_DECIDE;
 }

PETScLinearEquation:: ~PETScLinearEquation()
{
  VecDestroy(&b);
  VecDestroy(&x);
  MatDestroy(&A);
  if(lsolver) KSPDestroy(&lsolver);
  // if(prec) PCDestroy(&prec);

  if(global_x0)
    delete []  global_x0;
  if(global_x1)
    delete []  global_x1;
  if(global_buff)
    delete []  global_buff;

  PetscPrintf(PETSC_COMM_WORLD,"\n>>Number of Unknows: %d", m_size);
  PetscPrintf(PETSC_COMM_WORLD,"\n>>Elapsed time in linear solver: %fs", time_elapsed);
}

void PETScLinearEquation::Init(const int size, const int *sparse_index)
{
   m_size = size;

   if(sparse_index)
   {
      d_nz = sparse_index[0]; 
      o_nz = sparse_index[1]; 	
      nz = sparse_index[2]; 
      m_size_loc = sparse_index[3];          
   }   
    
   VectorCreate(m_size);   
   MatrixCreate(m_size, m_size);

   global_x0 = new PetscScalar[m_size];
   global_x1 = new PetscScalar[m_size];
   global_buff = new PetscScalar[m_size];

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
void PETScLinearEquation::Config(boost::property_tree::ptree &option)

{

   int maxits = 0;


   boost::optional<ptree> ptSolver = option.get_child("LinearSolver");
    if (!ptSolver)
        return;

    boost::optional<std::string> solver_type = ptSolver->get_optional<std::string>("solver_type");
    if (solver_type) {
      sol_type  = *solver_type;
    }
    boost::optional<std::string> precon_type = ptSolver->get_optional<std::string>("precon_type");
    if (precon_type) {
      pc_type = *precon_type;
    }
    boost::optional<double> error_tolerance = ptSolver->get_optional<double>("error_tolerance");
    if (error_tolerance) {
        ltolerance = *error_tolerance;
    }
    boost::optional<int> max_iteration_step = ptSolver->get_optional<int>("max_iteration_step");
    if (max_iteration_step) {
        maxits = *max_iteration_step;
    }

    const KSPType lsol_type = sol_type.c_str();
    const PCType prec_type = pc_type.c_str(); 


   KSPCreate(PETSC_COMM_WORLD,&lsolver);
   KSPSetOperators(lsolver, A, A,DIFFERENT_NONZERO_PATTERN);
   KSPSetType(lsolver,lsol_type);

   KSPGetPC(lsolver, &prec);
   PCSetType(prec, prec_type); //  PCJACOBI); //PCNONE);
   KSPSetTolerances(lsolver,ltolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxits);
   KSPSetFromOptions(lsolver);



}
//-----------------------------------------------------------------
void PETScLinearEquation::VectorCreate(PetscInt m)
{
  //PetscErrorCode ierr;  // returned value from PETSc functions 
  VecCreate(PETSC_COMM_WORLD, &b);
  ////VecCreateMPI(PETSC_COMM_WORLD,m_size_loc, m, &b);
  //VecSetSizes(b, m_size_loc, m);
  VecSetSizes(b, PETSC_DECIDE, m);
  VecSetFromOptions(b);
  VecDuplicate(b, &x);

  //VecGetOwnershipRange(b, &i_start,&i_end);
}


void PETScLinearEquation::MatrixCreate( PetscInt m, PetscInt n)
{

  MatCreate(PETSC_COMM_WORLD, &A);
  // TEST  MatSetSizes(A, m_size_loc, PETSC_DECIDE, m, n);

  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);
  //MatSetSizes(A, m_size_loc, PETSC_DECIDE, m,  n);

  MatSetFromOptions(A);
  MatSetType(A,MATMPIAIJ);
  MatMPIAIJSetPreallocation(A,d_nz,PETSC_NULL, o_nz,PETSC_NULL);
  //MatSeqAIJSetPreallocation(A,d_nz,PETSC_NULL);
  MatGetOwnershipRange(A,&i_start,&i_end);

  // std::cout<<"sub_a  "<<i_start<<";   sub_d "<<i_end<<std::endl;
}

void  PETScLinearEquation::getLocalRowColumnSizes(int *m, int *n)
{
  MatGetLocalSize(A, m, n);
}
void  PETScLinearEquation::getOwnerRange(int *start_r, int *end_r)
{
  *start_r = i_start;
  *end_r = i_end;
}

void PETScLinearEquation::Solver()
{
  
   //TEST
#ifdef TEST_MEM_PETSC
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif
 
  /* 
  //TEST
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.txt", &viewer);
  PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  PetscObjectSetName((PetscObject)x,"Solution");
  VecView(x, viewer);   
  */


   int its; 
   PetscLogDouble v1,v2;
   KSPConvergedReason reason;

   PetscGetTime(&v1);

   KSPSolve(lsolver, b, x);
  
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

   
#define aTEST_OUT
#ifdef TEST_OUT
  //TEST
   PetscViewer viewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x2.txt", &viewer);
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
   PetscObjectSetName((PetscObject)A,"Matrix");
   MatView(A, viewer);
   PetscObjectSetName((PetscObject)x,"Solution");
   VecView(x, viewer);
   PetscObjectSetName((PetscObject)b,"RHS");
   VecView(b, viewer);   
    VecDestroy(&b);
  VecDestroy(&x);
  MatDestroy(&A);
  if(lsolver) KSPDestroy(&lsolver);
  // if(prec) PCDestroy(&prec);
  if(global_x0)
    delete []  global_x0;
  if(global_x1)
    delete []  global_x1;
   PetscFinalize();
   exit(0);
#endif


#ifdef TEST_MEM_PETSC
  //TEST
   PetscMemoryGetCurrentUsage(&mem2);
   PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

  void PETScLinearEquation::AssembleRHS_PETSc()
{
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
}
void PETScLinearEquation::AssembleUnkowns_PETSc()
{
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
}
void PETScLinearEquation::AssembleMatrixPETSc(const MatAssemblyType type)
{
  MatAssemblyBegin(A, type);
  MatAssemblyEnd(A, type);
}


void PETScLinearEquation::finalAssemble()
{
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

}



void PETScLinearEquation::applyKnownSolutions(PetscInt ni,const PetscInt ix[], const PetscScalar y[])
{
    VecSetValues(x, ni, ix, y, INSERT_VALUES); 
    VecSetValues(b, ni, ix, y, INSERT_VALUES); 

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    
    zeroRows_in_Matrix(ni, ix);
    AssembleMatrixPETSc();
} 

void PETScLinearEquation::updateSolutions(PetscScalar *u0, PetscScalar *u1)
{

#ifdef TEST_MEM_PETSC
   //TEST
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif


  int i, j;
  PetscScalar *xp;
 
  int receivecount;
  PetscInt low,high,otherlow;
  MPI_Status status; 
  PetscInt count;
  int tag = 9999;
  VecGetOwnershipRange(x, &low, &high);
  VecGetLocalSize(x, &count);


  VecGetArray(x, &xp);
  for(i=0; i<count; i++)
    u1[i] = xp[i];


  //double *global_buff = new double[m_size];


  // Collect solution from processes.
  for(j=0; j<count; j++)
    global_buff[low+j] = u1[j];
  for(i=0;i<mpi_size;i++) 
  {
     if(i != rank)
     {

       MPI_Sendrecv( &count, 1, MPI_INT, i,tag, 
                     &receivecount,1,MPI_INT,i,tag, PETSC_COMM_WORLD ,&status);
       MPI_Sendrecv( &low, 1, MPI_INT, i,tag,
                 &otherlow,1,MPI_INT,i,tag,PETSC_COMM_WORLD,&status );
       MPI_Sendrecv( u1, count, MPI_DOUBLE, i,tag,
                     u0,receivecount,MPI_DOUBLE,i,tag, PETSC_COMM_WORLD,&status  );
       for(j=0;j<receivecount;j++)
         global_buff[otherlow+j] = u0[j];
     }
  }


  //MPI_Barrier(PETSC_COMM_WORLD);
  // Copy the collected solution to the array for the previous solution
  for(i=0;i<m_size;i++)
  {
    u1[i] = global_buff[i];
    u0[i] = global_buff[i];
  }
 

  //delete [] global_buff;

  VecRestoreArray(x, &xp);


  //TEST
#ifdef TEST_MEM_PETSC
   PetscMemoryGetCurrentUsage(&mem2);
   PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif

}

void PETScLinearEquation::mappingSolution()
{
  updateSolutions(global_x0, global_x1);
}


int PETScLinearEquation::getLocalSolution(PetscScalar *x_l)
{
  PetscInt count;
  VecGetLocalSize(x, &count);

  VecGetArray(x, &x_l);

  return count;
}


int PETScLinearEquation::getLocalRHS(PetscScalar *rhs_l)
{
  PetscInt count;
  VecGetLocalSize(b, &count);

  VecGetArray(b, &rhs_l);

  return count;
}

double *PETScLinearEquation::getGlobalSolution() const
{
  return global_x1;
}

/*!
  Get values of the specified elements from a global vector

  @param v_type - Indicator for vector: 0: x; 1: rhs
  @param ni 	- number of elements to get
  @param ix 	- indices where to get them from (in global 1d numbering) 
*/
void  PETScLinearEquation::getVecValues(const int v_type, PetscInt ni,const PetscInt ix[], 
				      PetscScalar y[]) const
{
  if(v_type == 0)
    VecGetValues(x, ni, ix, y);
  else 
    VecGetValues(b, ni, ix, y);
}

/*!
    Get norm of RHS
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i| 
    06.2012. WW
*/
PetscReal PETScLinearEquation::getVecNormRHS(NormType  nmtype)
{
  PetscReal norm = 0.;
  VecNorm(b, nmtype, &norm); 
  return norm; 
}
/*!
    Get norm of x
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i| 
    06.2012. WW
*/
PetscReal PETScLinearEquation::getVecNormX(NormType  nmtype)
{
  PetscReal norm = 0.;
  VecNorm(x, nmtype, &norm); 
  return norm; 
}


void  PETScLinearEquation::restoreLocalSolutionArray(PetscScalar *x_l)
{
   VecRestoreArray(x, &x_l);
}
void  PETScLinearEquation::restoreLocalRHSArray(PetscScalar *rhs_l)
{
   VecRestoreArray(b, &rhs_l);
}

void PETScLinearEquation::set_bVectorEntry(const int i, const double value )
{

  VecSetValues(b,1,&i,&value,INSERT_VALUES);
}
void PETScLinearEquation::set_xVectorEntry(const int i, const double value)
{

  VecSetValues(x,1,&i,&value,INSERT_VALUES);
}

void  PETScLinearEquation::setArrayValues(int arr_idx, PetscInt ni, const PetscInt ix[], 
                                       const PetscScalar y[],InsertMode iora) 
{
   if(arr_idx == 0)
     VecSetValues(x, ni, ix, y, iora); 
   else if(arr_idx == 1)
     VecSetValues(b, ni, ix, y, iora); 
}



void PETScLinearEquation::add_bVectorEntry(const int i, const double value,InsertMode mode )
{

  VecSetValue(b, i, value, mode);
}
void PETScLinearEquation::add_xVectorEntry(const int i, const double value, InsertMode mode)
{

  VecSetValue(x, i, value,mode);
}


void PETScLinearEquation::initializeMatVec( )
{

   VecSet(b, 0.0);
   VecSet(x, 0.0);
   MatZeroEntries(A);
} 



void PETScLinearEquation::addMatrixEntry(const int i, const int j, const double value)
{

  MatSetValue(A, i, j, value, ADD_VALUES);
}

void PETScLinearEquation::addMatrixEntries(const int m,const int idxm[], const int n, 
             const int idxn[],const PetscScalar v[])
{

  MatSetValues(A, m, idxm, n, idxn, v, ADD_VALUES);
}




void PETScLinearEquation::zeroRows_in_Matrix(const int nrows, const  PetscInt *rows)
{
  PetscScalar one = 1.0;
  // Each process indicates only rows it owns that are to be zeroed
  // MatSetOption(A, MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
  if(nrows>0)
    MatZeroRows (A, nrows, rows, one, PETSC_NULL, PETSC_NULL);
  else
    MatZeroRows(A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}

void PETScLinearEquation::Viewer(std::string file_name)
{
  PetscViewer viewer;
  std::string fname = file_name + "_eqs_dump.txt";
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
  PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);

  AssembleRHS_PETSc();
  AssembleUnkowns_PETSc();
  AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY );


  //PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
  PetscObjectSetName((PetscObject)A,"Stiffness_matrix");
  PetscObjectSetName((PetscObject)b,"RHS");
  PetscObjectSetName((PetscObject)x,"Solution");
  MatView(A,viewer);
  VecView(b, viewer);
  VecView(x, viewer);  

#define  nEXIT_TEST 
#ifdef EXIT_TEST 
  VecDestroy(&b);
  VecDestroy(&x);
  MatDestroy(&A);
  if(lsolver) KSPDestroy(&lsolver);
  // if(prec) PCDestroy(&prec);
  if(global_x0)
    delete []  global_x0;
  if(global_x1)
    delete []  global_x1;
   PetscFinalize();
   exit(0);
#endif 
 
}

} //end of namespace
 
