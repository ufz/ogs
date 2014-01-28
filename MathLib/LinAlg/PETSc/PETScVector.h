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

        PETScVector() = default;
        /*!
             \brief  Constructor
             \param size  the global size of the vector
         */
        explicit PETScVector(const PetscInt size)
        {
            _size = size;
            create(_size);
        }

        /// copy constructor
        PETScVector(const PETScVector &existing_vec);

        /// Duplicate a vector. The value of existing_vec is not copied.
        void Duplicate(const PETScVector &existing_vec);

        ~PETScVector()
        {
            VecDestroy(&_v);
        }

        /// Perform MPI collection of assembled entries in buffer
        void finalizeAssembly()
        {
            VecAssemblyBegin(_v);
            VecAssemblyEnd(_v);
        }

        /// Get the global size the vector
        PetscInt  size() const
        {
            return _size;
        }

        /// Get the number of entries in the same rank
        PetscInt  getLocalSize() const
        {
            return _size_loc;
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

        /*!
          Get norm of vector
          \param nmtype   norm type
                          sum_abs_entries denotes \f$\sum_i |x_i|\f$
                          euclidean denotes \f$\sqrt(\sum_i (x_i)^2)\f$
                          max_abs_entry denotes \f$\mathrm{max}_i |x_i|\f$

           \var VectorNormType  an enum type has values as
                     SUM_ABS_ENTRIES = 0,
                     EUCLIDEAN = 1,
                     MAX_ABS_ENTRY = 2
        */
        PetscScalar getNorm(const MathLib::VecNormType nmtype = MathLib::VecNormType::NORM2) const;

        /*!
           Insert a single entry with value.

           \param i  entry index
           \param value  entry value

         */
        void set(const PetscInt i, const PetscScalar value)
        {
            VecSetValue(_v, i, value, INSERT_VALUES);
        }

        /*!
           Add a value to an entry.

           \param i  number of the entry
           \param value value.
         */
        void add(const PetscInt i, const PetscScalar value)
        {
            VecSetValue(_v, i, value,  ADD_VALUES);
        }

        /*!
           Add values to several entries
           \param e_idxs  indicies of entries to be added
                          Note: size_t cannot be the type of e_idxs template argument
           \param sub_vec entries to be added
        */
        template<class T_SUBVEC>  void add(const std::vector<PetscInt> &e_idxs,
                                           const T_SUBVEC &sub_vec)
        {
            VecSetValues(_v, e_idxs.size(), &e_idxs[0], &sub_vec[0], ADD_VALUES);
        }

        /*!
           Add values to several entries
           \param e_idxs  indicies of entries to be added.
                          Note: size_t cannot be the type of e_idxs template argument
           \param sub_vec entries to be added
        */
        template<class T_SUBVEC>  void set(const std::vector<PetscInt> &e_idxs,
                                           const T_SUBVEC &sub_vec)
        {
            VecSetValues(_v, e_idxs.size(), &e_idxs[0], &sub_vec[0], INSERT_VALUES);
        }

        /*!
           get several entries
           \param e_idxs  indicies of entries to be gotten.
                          Note: size_t cannot be the type of e_idxs template argument
           \param sub_vec values of entries
        */
        template<class T_SUBVEC>  void get(const std::vector<PetscInt> &e_idxs,
                                           T_SUBVEC &sub_vec)
        {
            VecGetValues(_v, e_idxs.size(), &e_idxs[0], &sub_vec[0]);
        }

        /*!
           Get local vector, i.e. entries in the same rank

           \param loc_vec  pinter to array where stores the local vector
         */
        void getLocalVector(PetscScalar *loc_vec) const
        {
            PetscInt count;
            VecGetLocalSize(_v, &count);
            VecGetArray(_v, &loc_vec);
        }

        /*!
           Get global vector

           \param u0  array to store the global vector
           \param u1  array to store the global vector too, and it is also used as
                      a buffer array in this function.
         */
        void getGlobalEntries(PetscScalar u0[], PetscScalar u1[]);

        /// Get an entry value. This is an expensive operation,
        /// and it only get local value. Use it for only test purpose
        PetscScalar get(const  PetscInt idx) const;

        /// Get PETsc vector. Use it only for test purpose
        PETSc_Vec &getData()
        {
            return _v;
        }

        /// Initialize  the vector with a constant value
        void operator = (const PetscScalar val)
        {
            VecSet(_v, val);
        }
        /// Overloaded operator: assign
        void operator = (const PETScVector &v_in)
        {
            VecCopy(_v, v_in._v);
        }

        ///  Overloaded operator: add
        void operator += (const PETScVector& v_in)
        {
            VecAXPY(_v, 1.0, v_in._v);
        }

        ///  Overloaded operator: subtract
        void operator -= (const PETScVector& v_in)
        {
            VecAXPY(_v, -1.0, v_in._v);
        }

        /*! View the global vector for test purpose. Do not use it for output a big vector.
            \param file_name  file name for output
            \vw_format  file format listed as:
        PETSC_VIEWER_DEFAULT 	- default format
        PETSC_VIEWER_ASCII_MATLAB 	- MATLAB format
        PETSC_VIEWER_ASCII_DENSE 	- print matrix as dense
        PETSC_VIEWER_ASCII_IMPL 	- implementation-specific format (which is in many cases the same as the default)
        PETSC_VIEWER_ASCII_INFO 	- basic information about object
        PETSC_VIEWER_ASCII_INFO_DETAIL 	- more detailed info about object
        PETSC_VIEWER_ASCII_COMMON 	- identical output format for all objects of a particular type
        PETSC_VIEWER_ASCII_INDEX 	- (for vectors) prints the vector element number next to each vector entry
        PETSC_VIEWER_ASCII_SYMMODU 	- print parallel vectors without indicating the processor ranges
        PETSC_VIEWER_ASCII_VTK 	- outputs the object to a VTK file
        PETSC_VIEWER_NATIVE 	- store the object to the binary file in its native format (for example, dense matrices are stored as dense), DMDA vectors are dumped directly to the file instead of being first put in the natural ordering
        PETSC_VIEWER_DRAW_BASIC 	- views the vector with a simple 1d plot
        PETSC_VIEWER_DRAW_LG 	- views the vector with a line graph
        PETSC_VIEWER_DRAW_CONTOUR 	- views the vector with a contour plot

        */
        void viewer(const std::string &file_name,
                    const PetscViewerFormat vw_format = PETSC_VIEWER_ASCII_MATLAB );

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

        /// Wrap the PETSc function VecRestoreArray to restore a vector after VecGetArray is called.
        void restoreLocalVector(PetscScalar *loc_vec)
        {
            VecRestoreArray(_v, &loc_vec);
        }

        /*!
              \brief  collect local vectors
              \param  local_array  local array
              \param  global_array global array
        */
        void collectLocalVectors(PetscScalar local_array[],
                                 PetscScalar global_array[]);
};

} // end namespace
#endif


