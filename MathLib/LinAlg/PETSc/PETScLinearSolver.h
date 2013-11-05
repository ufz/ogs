/*!
   \file  PETScLinearEquation.h
   \brief Declaration of class PETScLinearSolve, which defines a solver object
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
#ifndef PETScLINEAR_SOLVER_H_
#define PETScLINEAR_SOLVER_H_

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

#include "petscmat.h"
#include "petscksp.h"



namespace MathLib
{

class PETScMatrix;
class PETScVector;

class PETScLinearSolver
{
   public:
      explicit PETScLinearSolver(PETScMatrix &stiffness_matrix, PETScVector &rhs, PETScVector &unknowns);

      ~PETScLinearSolver();

      /// Set the solver properties
      void Config(boost::property_tree::ptree &option);

      void allocateMemory4TemoraryArrays(const PetscInt size);
      void releaseMemory4TemoraryArrays();


      void Solver();
      void Solver(PETScVector &rhs, PETScVector &unknowns);

      void finalAssemble();

      void addMatrixEntries(const PetscInt m,const PetscInt idxm[], const PetscInt n,
                            const PetscInt idxn[], const PetscScalar v[]);


      void applyKnownSolutions(PetscInt ni,const PetscInt ix[], const PetscScalar y[]);

      void mappingSolution();

      PetscScalar *getGlobalSolution() const;

      PetscReal getNormRHS(NormType  nmtype);
      PetscReal getNormUnknowns(NormType  nmtype);

      /// Get the number of global unknows
      int Size() const
      {
         return m_size;
      }

      //  PETScMatrix &getMatrix() { return A; }

      void initializeMatVec();


      void set_rank_size(const int mrank, const int msize)
      {
         mpi_size = msize;
         rank = mrank;
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
      PETScMatrix  &A;
      PETScVector  &b;
      PETScVector  &x;
      KSP lsolver;
      PC prec;

      /// Temporary arrays to store the global solutions
      PetscScalar *global_x0 = nullptr;
      PetscScalar  *global_x1 = nullptr;


      /// Slover and preconditioner names, only for log
      std::string sol_type;
      std::string pc_type;

      PetscLogDouble time_elapsed;

      /// Dimension of the unknows
      PetscInt m_size;
      /// Dimention of the local unknowns
      PetscInt m_size_loc;
      float ltolerance;

      /// Rank size
      int mpi_size;
      /// Rank
      int rank;

      /// Indictor whether A, b and x are all set
      //bool full_fill_A_b_x;
};

} // end namespace
#endif

