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
#include "LinAlg/VectorNorms.h"

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
        explicit PETScVector(const PetscInt size);

        /// Copy constructor. The value of existing_vec is not copied.
        PETScVector(const PETScVector &existing_vec);
        ~PETScVector();

        /// Perform MPI collection of assembled entries in buffer
        void finalizeAssembly();

        /// Get PETsc vector. Use it only for test purpose
        PETSc_Vec &getData()
        {
            return _v;
        }

        /*!
           Get local vector, i.e. entries in the same rank

           \param loc_vec  pinter to array where stores the local vector
           \return  the size of the local vector
         */
        PetscInt getLocalVector(PetscScalar *loc_vec) const;

        /// Wrap the PETSc function VecRestoreArray to restore a vector after VecGetArray is called.
        void restoreLocalVector(PetscScalar *loc_vec);

        /*!
         Get values of the specified elements from a global vector

         \param ni  number of elements to get
         \param ix  indices where to get them from (in global 1d numbering)
         \param y   values of the got entries
        */
        void getEntries(PetscInt ni, const PetscInt ix[],
                        PetscScalar y[]) const;
        /*!
           Get global vector

           \param u0  array to store the global vector
           \param u1  array to store the global vector too, and it is also used as
                      a buffer array in this function.
         */
        void getGlobalEntries(PetscScalar u0[], PetscScalar u1[]) const;

        /*!
          Get norm of vector
          \param nmtype   norm type
                          sum_abs_entries denotes \f$\sum_i |x_i|\f$
                          euclidean denotes \f$\sqrt(\sum_i (x_i)^2)\f$
                          max_abs_entry denotes \f$\mathrm{max}_i |x_i|\f$
        */
        PetscReal getNorm(const VectorNormType nmtype = EUCLIDEAN) const;

        /// Get the size the vector
        PetscInt  size() const
        {
            return _size;
        }

        /// Get the number of entries in the same rank
        PetscInt  getLocalSize() const
        {
            return _size_loc;
        }

        /*!
           Insert a single entry with value.

           \param i  entry index
           \param value  entry value

         */
        void set(const PetscInt i, const PetscScalar value);

        /*!
           Add a value to an entry.

           \param i  number of the entry
           \param value value.
         */
        void add(const PetscInt i, const PetscScalar value);

        /*!
           Insert multi-entries with value.

           \param ni number of entries
           \param ix indicies of entries
           \param y  values of entries
         */
        void setValues(PetscInt ni,const PetscInt ix[],
                       const PetscScalar y[]);

        /*!
           Add values to multi-entries.

           \param ni number of entries
           \param ix indicies of entries
           \param y  values of entries
         */
        void addValues(PetscInt ni,const PetscInt ix[],
                       const PetscScalar y[]);

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
                this->add(e_idxs[i], sub_vec[i]);
            }
        }

        /// Set all entries to zero value
        void setZero();

        /// Get an entry value. This is an expensive operation,
        /// and it only get local value. Use it for only test purpose
        double get(const  PetscInt idx) const;

        /// Initialize  the vector with a constant value
        void operator = (const PetscScalar val);

        /// Overloaded operator: assign
        void operator = (const PETScVector &v_in);

        ///  Overloaded operator: add
        void operator += (const PETScVector& v_in);

        ///  Overloaded operator: subtract
        void operator -= (const PETScVector& v_in);

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

        /// View the global vector for test purpose. Do not use it for output a big vector.
        void viewer(const std::string &file_name);

    private:
        PETSc_Vec _v;

        /// Starting index in a rank
        PetscInt _start_rank;
        /// Ending index in a rank
        PetscInt _end_rank;

        /// Size of the vector
        PetscInt _size;
        /// Size of local entries
        PetscInt _size_loc;

        /// Create vector
        void create(PetscInt vec_size);
};

} // end namespace
#endif


