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
      explicit PETScVector(const PetscInt size);

      /// Copy constructor. The value of existing_vec is not copied.
      PETScVector(const PETScVector &existing_vec);
      ~PETScVector();


      /*!
          \brief Set the number of unknowns, and create vector. Iddentical to PETScVector(const PetscInt size);

          \param vec_size,  the number of unknowns.

      */
      void Init(const PetscInt vec_size);

      /// Perform MPI collection of assembled entries in buffer
      void finalAssemble();

      PETSc_Vec &getData()
      {
         return v;
      }


      int getLocalVector(PetscScalar *loc_vec);


      /*!
       Get values of the specified elements from a global vector

       \param v_type - Indicator for vector: 0: x; 1: rhs
       \param ni 	- number of elements to get
       \param ix 	- indices where to get them from (in global 1d numbering)
      */     void getEntries(PetscInt ni, const PetscInt ix[],
                             PetscScalar y[]) const;


      void getGlobalEntries(PetscScalar *u0, PetscScalar *u1);


      /*!
        Get norm of vector
        \param nmtype  - norm type
                         NORM_1 denotes \f$\sum_i |x_i|\f$
                         NORM_2 denotes \f$\sqrt(\sum_i (x_i)^2)\f$
                         NORM_INFINITY denotes \f$\mathrm{max}_i |x_i|\f$
      */
      PetscReal getNorm(NormType nmtype);
      PetscReal getNorm()
      {
         return  getNorm(NORM_2);
      }

      void restoreLocalVector(PetscScalar *loc_vec);


      void getOwnerRange(int *start_r, int *end_r);

      /// Get the number of global unknows
      int size() const
      {
         return _size;
      }

      void set(const PetscInt i, const double value);

      void setValues( PetscInt ni,const PetscInt ix[],
                      const PetscScalar y[], InsertMode iora = ADD_VALUES);

      void add(const PetscInt i, const double value, InsertMode mode = ADD_VALUES);


      /*!
         Add values to several entries
         \e_idxs   - Indicies of entries to be added
         \sub_vec  - entries to be added
      */

      template<class T_SUBVEC>  void add(const std::vector<std::size_t> &e_idxs,
                                         const T_SUBVEC &sub_vec)
      {
         for (std::size_t i=0; i<e_idxs.size(); ++i)
         {
            add(e_idxs[i], sub_vec[i], ADD_VALUES);
         }
      }

      void setZero();


      /// Get an entry value. This is an expensive operation,
      /// and it only get local value. Use it for only test purpose
      double get(const  PetscInt idx) const;


      /// Initialize  the vector with a constant value
      void operator = (const double val);

      // Overloaded operators:
      /// Overloaded operator: asign
      void operator = (const PETScVector &v_in);

      /// Overloaded operator: assignment
      PETScVector& operator = (PETScVector &v_in);

      ///  Overloaded operator: add
      void operator += (const PETScVector& v_in);

      ///  Overloaded operator: subtract
      void operator -= (const PETScVector& v_in);



      void set_rank_size(const int myrank, const int ranksize)
      {
         _size_rank = ranksize;
         rank = myrank;
      }


      PetscInt getRangeBegin() const
      {
         return _start_rank;
      }
      PetscInt getRangeEnd() const
      {
         return _end_rank;
      }

      PetscInt getMPI_Size() const
      {
         return _size_rank;
      }
      PetscInt getMPI_Rank() const
      {
         return rank;
      }

      void Viewer(std::string file_name);


   private:
      PETSc_Vec v;

      /// Starting index in a rank
      PetscInt _start_rank;
      /// Ending index in a rank
      PetscInt _end_rank;

      /// Dimension of the unknows
      PetscInt _size;
      /// Dimention of the local unknowns
      PetscInt _size_loc;

      /// Rank size
      int _size_rank;
      /// Rank
      int rank;

      /// Create vector
      void Create(PetscInt vec_size);

};

} // end namespace
#endif

