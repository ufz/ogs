/*!
   \file  PETScLinearEquation.h
   \brief Declaration of class PETScLinearEquation, which provides interfaces to
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
#ifndef PETScLINEAR_EQUATION_H_
#define PETScLINEAR_EQUATION_H_

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>


#include "petscmat.h"
#include "petscksp.h"

typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;

namespace MathLib
{

class PETScLinearEquation
{
   public:
      PETScLinearEquation();
      ~PETScLinearEquation();

      /// Set the solver properties
      void Config(boost::property_tree::ptree &option);

      /*!
          \brief Set the number of unknowns, sparsity of matrix, and create matrix and vector

          @param size,  the number of unknowns.
          @param sparse_index,  the sparsity index.

       */
      void Init(const int size, const int *sparse_index = NULL);

      void Solver();
      void AssembleRHS_PETSc();
      void AssembleUnkowns_PETSc();
      void AssembleMatrixPETSc(const MatAssemblyType type = MAT_FINAL_ASSEMBLY); //MAT_FLUSH_ASSEMBLY
      void finalAssemble();
      void applyKnownSolutions(PetscInt ni,const PetscInt ix[], const PetscScalar y[]);

      void updateSolutions(PetscScalar *u0, PetscScalar *u1);
      void mappingSolution();
      int getLocalSolution(PetscScalar *x_l);
      int getLocalRHS(PetscScalar *rhs_l);
      double *getGlobalSolution() const;
      void  getVecValues(const int v_type, PetscInt ni, const PetscInt ix[],
                         PetscScalar y[]) const;
      PetscReal getVecNormRHS(NormType  nmtype= NORM_2);
      PetscReal getVecNormX(NormType  nmtype= NORM_2);

      void restoreLocalSolutionArray(PetscScalar *x_l);
      void restoreLocalRHSArray(PetscScalar *rhs_l);
      void getLocalRowColumnSizes(int *m, int *n);
      void getOwnerRange(int *start_r, int *end_r);

      /// Get the number of global unknows
      int Size() const
      {
         return m_size;
      }

      void set_xVectorEntry(const int i, const double value);
      void set_bVectorEntry(const int i, const double value);
      void setArrayValues(int arr_idx, PetscInt ni,const PetscInt ix[],
                          const PetscScalar y[], InsertMode iora = ADD_VALUES);

      void add_xVectorEntry(const int i, const double value, InsertMode mode);
      void add_bVectorEntry(const int i, const double value, InsertMode mode);
      void addMatrixEntry(const int i, const int j, const double value);
      void addMatrixEntries(const int m,const int idxm[], const int n,
                            const int idxn[],const PetscScalar v[]);

      void initializeMatVec();

      void zeroRows_in_Matrix(const int nrow, const  PetscInt *rows);
      void zeroMatrix()
      {
         MatZeroEntries(A);
      }


      void set_rank_size(const int mrank, const int msize)
      {
         mpi_size = msize;
         rank = mrank;
      }


      PetscInt getStartRow() const
      {
         return i_start;
      }
      PetscInt getEndRow() const
      {
         return i_end;
      }

      PetscInt getMPI_Size() const
      {
         return mpi_size;
      }
      PetscInt getMPI_Rank() const
      {
         return rank;
      }

      void Viewer(std::string file_name);

   private:
      PETSc_Mat  A;
      PETSc_Vec b;
      PETSc_Vec x;
      KSP lsolver;
      PC prec;
      PetscInt i_start;
      PetscInt i_end;

      PetscScalar *global_x0;
      PetscScalar *global_x1;
      PetscScalar *global_buff;

      // Slover and preconditioner names, only for log
      std::string sol_type;
      std::string pc_type;

      PetscLogDouble time_elapsed;

      /// Dimension of the unknows
      PetscInt m_size;
      /// Dimention of the local unknowns
      PetscInt m_size_loc;
      float ltolerance;
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

      /// Create vector
      void VectorCreate(PetscInt m);
      /// Creat matrix
      void MatrixCreate(PetscInt m, PetscInt n);

};

// extern std::vector<PETScLinearEquation*> EQS_Vector;
} // end namespace
#endif

