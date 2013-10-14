/*!
   \file  PETScMatrix.h
   \brief Declaration of class PETScMatrix, which provides an interface to
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
#ifndef PETScMATRIX_H_
#define PETScMATRIX_H_
 
#include <string>
#include <vector>


#include "petscmat.h"
#include "petscksp.h"

typedef Mat PETSc_Mat;

namespace MathLib
{
 
class PETScMatrix
{
public:    
    PETScMatrix();
    ~PETScMatrix();

    void Init(const PetscInt size, const int *sparse_index = NULL);
  
  
    void getLocalRowColumnSizes(int *m, int *n); 
    void getOwnerRange(int *start_r, int *end_r); 

    /// Get the number of global unknows
    int Size() const {return m_size;}


    /// Allocate mmemory for matrix
    void Create(const  PetscInt m, const PetscInt n);

    void zeroRows_in_Matrix(const PetscInt nrow, const  PetscInt *rows);
    void zeroMatrix()
    {
       MatZeroEntries(A);
    }

    void addEntry(const PetscInt i, const PetscInt j, const PetscScalar value);

    /*!
      \brief Add a submatrix to this
    
      @param m,     number of rows of the submatrix.
      @param idxm,  global row indicies of the submatrix
      @param n,     number of columns of the submatrix.
      @param idxn,  global culumn indicies of the submatrix.
      @param v,     valuse.
   */

    void addEntries(const PetscInt m, const PetscInt idxm[], 
                    const PetscInt n, const PetscInt idxn[], 
                    const PetscScalar v[]);

    /// Perform MPI collection of assembled entries in buffer
    void finalAssemble(const MatAssemblyType type = MAT_FINAL_ASSEMBLY); //MAT_FLUSH_ASSEMBLY
 


    void set_rank_size(const int mrank, const int msize)
    {
      mpi_size = msize;
      rank = mrank;  
    } 
 
    PETSc_Mat &getData() 
    {
      return A;  
    } 
    

    PetscInt getStartRow() const {return i_start;} 
    PetscInt getEndRow() const {return i_end;} 

    PetscInt getMPI_Size() const {return mpi_size;} 
    PetscInt getMPI_Rank() const {return rank;} 
   
    void Viewer(std::string file_name);

  private:
    PETSc_Mat  A;
    /// Starting index in a rank
    PetscInt i_start;
    /// Ending index in a rank
    PetscInt i_end;

    /// Dimension of the unknows
    PetscInt m_size;
    /// Dimention of the local unknowns
    PetscInt m_size_loc;

    /// Number of nonzeros per row in DIAGONAL portion of 
    /// local submatrix (same value is used for all local rows) 
    PetscInt d_nz; 
    /// Number of nonzeros per row in the OFF-DIAGONAL portion of 
    /// local submatrix (same value is used for all local rows). 
    PetscInt o_nz; 	
    /// Number of nonzeros per row (same for all rows) 
    PetscInt nz;

    /// Rank size
    int mpi_size;
    /// Rank
    int rank;

    /// Creat matrix
    void MatrixCreate(PetscInt m, PetscInt n);

};

} // end namespace 
#endif

