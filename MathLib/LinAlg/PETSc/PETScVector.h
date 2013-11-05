/*!
   \file  PETScVector.h
   \brief Declaration of class PETScVector, which provides an interface to
          PETSc vector routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013


  \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/
#ifndef PETScVECTOR_H_
#define PETScVECTOR_H_

#include <string>
#include <vector>

#include "petscvec.h"

typedef Vec PETSc_Vec;

namespace MathLib
{

class PETScVector
{
   public:
      PETScVector();
      PETScVector(const PetscInt size);
      PETScVector(PETScVector &existing_vec);
      ~PETScVector();


      /*!
          \brief Set the number of unknowns, and create vector. Iddentical to PETScVector(const PetscInt size);

          @param size,  the number of unknowns.

       */
      void Init(const PetscInt size);

      /// Perform MPI collection of assembled entries in buffer
      void finalAssemble();

      PETSc_Vec &getData()
      {
         return v;
      }


      int getLocalVector(PetscScalar *loc_vec);
      void getEntries(PetscInt ni, const PetscInt ix[],
                      PetscScalar y[]) const;


      void getGlobalEntries(PetscScalar *u0, PetscScalar *u1);


      /*!
        Get norm of vector
        @param nmtype  - norm type
                         NORM_1 denotes \f$\sum_i |x_i|\f$
                         NORM_2 denotes \f$\sqrt(\sum_i (x_i)^2)\f$
                         NORM_INFINITY denotes \f$\mathrm{max}_i |x_i|\f$
         06.2012.
      */
      PetscReal getNorm(NormType  nmtype= NORM_2);

      void restoreLocalVector(PetscScalar *loc_vec);


      void getOwnerRange(int *start_r, int *end_r);

      /// Get the number of global unknows
      int size() const
      {
         return m_size;
      }

      void set(const PetscInt i, const double value);

      void setValues( PetscInt ni,const PetscInt ix[],
                      const PetscScalar y[], InsertMode iora = ADD_VALUES);

      void addValue(const PetscInt i, const double value, InsertMode mode);

      /// Add values to several entries
      template<class T_SUBVEC>
      void add(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec)
      {
         for (std::size_t i=0; i<pos.size(); ++i)
         {
            addValue(pos[i], sub_vec[i], ADD_VALUES);
         }
      }

      void nullize();



      void set_rank_size(const int mrank, const int msize)
      {
         mpi_size = msize;
         rank = mrank;
      }


      PetscInt getRangeStart() const
      {
         return i_start;
      }
      PetscInt getRangeEnd() const
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
      PETSc_Vec v;

      /// Starting index in a rank
      PetscInt i_start;
      /// Ending index in a rank
      PetscInt i_end;

      /// Dimension of the unknows
      PetscInt m_size;
      /// Dimention of the local unknowns
      PetscInt m_size_loc;

      /// Rank size
      int mpi_size;
      /// Rank
      int rank;

      /// Create vector
      void Create(PetscInt m);

};

} // end namespace
#endif

