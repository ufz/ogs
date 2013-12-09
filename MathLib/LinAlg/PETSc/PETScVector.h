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
#ifndef PETSCVECTOR_H_
#define PETSCVECTOR_H_

#include <string>
#include <vector>

#include "petscvec.h"

typedef Vec PETSc_Vec;

namespace MathLib
{
/*!
   \class PETScVector

   \brief Wrapper class for PETSc vector
*/
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

          \param vec_size  the number of unknowns.

      */
      void init(const PetscInt vec_size);

      /// Perform MPI collection of assembled entries in buffer
      void finalizeAssembly();

      /// Get PETsc vector. Only used it for test purpose
      PETSc_Vec &getData()
      {
         return v;
      }

      /*!
         Get local vector, i.e. entries in the same rank

         \param loc_vec  array to store the local vector
         \return  the size of the local vector
       */
      PetscInt getLocalVector(PetscScalar loc_vec[]) const;

      /*!
       Get values of the specified elements from a global vector

       \param v_type  indicator for vector: 0: x; 1: rhs
       \param ni  number of elements to get
       \param ix  indices where to get them from (in global 1d numbering)
      */
      void getEntries(PetscInt ni, const PetscInt ix[],
                      PetscScalar y[]) const;

      /*!
         Get global vector

         \param u0  array to store the global vector
         \param u1  array to store the global vector too, and it is also used as
                    a buffer array in this function.
       */
      void getGlobalEntries(PetscScalar u0[], PetscScalar u1[]);

      /*!
        Get norm of vector
        \param nmtype   norm type
                         NORM_1 denotes \f$\sum_i |x_i|\f$
                         NORM_2 denotes \f$\sqrt(\sum_i (x_i)^2)\f$
                         NORM_INFINITY denotes \f$\mathrm{max}_i |x_i|\f$
      */
      PetscReal getNorm(NormType nmtype) const;

      /// Special case of getNorm(NormType nmtype), only for the consistency of the ogs6 vector template.
      PetscReal getNorm() const
      {
         return  getNorm(NORM_2);
      }

      /// Wrap the PETSc function VecRestoreArray to restore a vector after VecGetArray is called.
      void restoreLocalVector(PetscScalar loc_vec[]);

      /// Get the number of global unknows
      int size() const
      {
         return _size;
      }

      /*!
         Insert a single entry with value.

         \param i  entry index
         \param value  entry value

       */
      void set(const PetscInt i, const PetscScalar value);

      /*!
         Insert multi-entries with value, or add values to multi-entries.

         \param ni number of entries
         \param ix indicies of entries
         \param y  values of entries
         \param iora  indicator of action: add or insert

       */
      void setValues(PetscInt ni,const PetscInt ix[],
                     const PetscScalar y[], InsertMode mode = ADD_VALUES);

      /*!
         Insert  an entry with value, or add a value to an entry.

         \param i  number of the entry
         \param value value.
         \param mode indicator of action: add or insert

       */
      void add(const PetscInt i, const PetscScalar value, InsertMode mode = ADD_VALUES);

      /*!
         Add values to several entries
         \e_idxs  indicies of entries to be added
         \sub_vec entries to be added
      */

      template<class T_SUBVEC>  void add(const std::vector<std::size_t> &e_idxs,
                                         const T_SUBVEC &sub_vec)
      {
         for (std::size_t i=0; i<e_idxs.size(); ++i)
         {
            add(e_idxs[i], sub_vec[i], ADD_VALUES);
         }
      }

      /// Set all entries to zero value
      void setZero();

      /// Get an entry value. This is an expensive operation,
      /// and it only get local value. Use it for only test purpose
      double get(const  PetscInt idx) const;

      /// Initialize  the vector with a constant value
      void operator = (const PetscScalar val);

      /// Overloaded operator: asign
      void operator = (const PETScVector &v_in);

      /// Overloaded operator: assignment
      PETScVector& operator = (PETScVector &v_in);

      ///  Overloaded operator: add
      void operator += (const PETScVector& v_in);

      ///  Overloaded operator: subtract
      void operator -= (const PETScVector& v_in);

      /// Get rank index and the number of processors
      void set_rank_size(const int myrank, const int ranksize)
      {
         _size_rank = ranksize;
         _rank = myrank;
      }

      /// Get the start index of the local vector
      PetscInt getRangeBegin() const
      {
         return _start_rank;
      }

      /// Get the end index of the local vector
      PetscInt getRangeEnd() const
      {
         return _end_rank;
      }

      /// Get the maximum rank
      PetscInt getSizeMPI() const
      {
         return _size_rank;
      }

      /// Get the current rank index
      PetscInt getRankMPI() const
      {
         return _rank;
      }

      /// View the global vector for test purpose. Do not use it for output a big vector.
      void Viewer(const std::string &file_name);

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
      int _rank;

      /// Create vector
      void create(PetscInt vec_size);
};

} // end namespace
#endif

