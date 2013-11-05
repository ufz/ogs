/*!
   \file  PETScMatrix.cpp
   \brief Definition of member functions of class PETScMatrix, which provides an interface to
          PETSc matrix routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013


  \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/


#include "PETScMatrix.h"

#include<iostream>


namespace MathLib
{


PETScMatrix :: PETScMatrix ()
{
   d_nz = 10;
   o_nz = 10;
   nz = 10;
   m_size_loc = PETSC_DECIDE;
}

PETScMatrix:: ~PETScMatrix()
{
   MatDestroy(&A);
}

void PETScMatrix::Init(const PetscInt size, const int *sparse_index)
{
   m_size = size;

   if(sparse_index)
   {
      d_nz = sparse_index[0];
      o_nz = sparse_index[1];
      nz = sparse_index[2];
      m_size_loc = sparse_index[3];
   }

   Create(m_size, m_size);
}


void PETScMatrix::finalAssemble(const MatAssemblyType type)
{
   MatAssemblyBegin(A, type);
   MatAssemblyEnd(A, type);
}


void PETScMatrix::Create(const  PetscInt m, const PetscInt n)
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

   //  std::cout<<"sub_a  "<<i_start<<";   sub_d "<<i_end<<std::endl;
}

void  PETScMatrix::getLocalRowColumnSizes(int *m, int *n)
{
   MatGetLocalSize(A, m, n);
}
void  PETScMatrix::getOwnerRange(int *start_r, int *end_r)
{
   *start_r = i_start;
   *end_r = i_end;
}

void PETScMatrix::setValue(const PetscInt i, const PetscInt j, const PetscScalar value)
{

   MatSetValue(A, i, j, value, INSERT_VALUES);
}


void PETScMatrix::add(const PetscInt i, const PetscInt j, const PetscScalar value)
{

   MatSetValue(A, i, j, value, ADD_VALUES);
}

void PETScMatrix::addEntries(const PetscInt m,const int idxm[], const PetscInt n,
                             const PetscInt idxn[],const PetscScalar v[])
{

   MatSetValues(A, m, idxm, n, idxn, v, ADD_VALUES);
}



void PETScMatrix::zeroRows_in_Matrix(const PetscInt nrows, const  PetscInt *rows)
{
   PetscScalar one = 1.0;
   // Each process indicates only rows it owns that are to be zeroed
   // MatSetOption(A, MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
   if(nrows>0)
      MatZeroRows (A, nrows, rows, one, PETSC_NULL, PETSC_NULL);
   else
      MatZeroRows(A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}



void PETScMatrix::Viewer(std::string file_name)
{
   PetscViewer viewer;
   std::string fname = file_name + "_stifness_matrix.txt";
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);

   finalAssemble(MAT_FINAL_ASSEMBLY );


   //PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
   PetscObjectSetName((PetscObject)A,"Stiffness_matrix");
   MatView(A,viewer);

#define  nEXIT_TEST
#ifdef EXIT_TEST
   MatDestroy(&A);
   PetscFinalize();
   exit(0);
#endif

}



bool finalizeMatrixAssembly(PETScMatrix &mat)
{
   mat.finalAssemble();
   return true;
}

} //end of namespace

