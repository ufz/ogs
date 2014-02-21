/**
 * \file  PETScVector.h
 * \brief Declaration of class PETScVector, which provides an interface to
 *        PETSc vector routines.
 *
 * \author Wenqing Wang
 * \date Nov 2011 - Sep 2013
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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

        /*!
            \brief Constructor
            \param vec_size       The size of the vector, either global or local
            \param is_global_size The flag of the global size, the default is true
        */
        PETScVector(const PetscInt vec_size, const bool is_global_size = true);

        /*!
             \brief Copy constructor
             \param existing_vec The vector to be copied
             \param deep_copy    The flag for a deep copy, which means to copy the values as well,
                                 the default is true
        */
        PETScVector(const PETScVector &existing_vec, const bool deep_copy = true);

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

        /// Get the global size of the vector
        PetscInt size() const
        {
            return _size;
        }

        /// Get the number of entries in the same rank
        PetscInt getLocalSize() const
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
          \param nmtype Norm type, default Euclidean norm
        */
        PetscScalar getNorm(const MathLib::VecNormType nmtype = MathLib::VecNormType::NORM2) const;

        /*!
           Insert a single entry with value.

           \param i     Entry index
           \param value Entry value

        */
        void set(const PetscInt i, const PetscScalar value)
        {
            VecSetValue(_v, i, value, INSERT_VALUES);
        }

        /*!
           Add a value to an entry.

           \param i     Number of the entry
           \param value Value.
        */
        void add(const PetscInt i, const PetscScalar value)
        {
            VecSetValue(_v, i, value,  ADD_VALUES);
        }

        /*!
           Add values to several entries
           \param e_idxs  Indicies of entries to be added
                          Note: size_t cannot be the type of e_idxs template argument
           \param sub_vec Entries to be added
        */
        template<class T_SUBVEC> void add(const std::vector<PetscInt> &e_idxs,
                                          const T_SUBVEC &sub_vec)
        {
            VecSetValues(_v, e_idxs.size(), &e_idxs[0], &sub_vec[0], ADD_VALUES);
        }

        /*!
           Add values to several entries
           \param e_idxs  Indicies of entries to be added.
                          Note: size_t cannot be the type of e_idxs template argument
           \param sub_vec Entries to be added
        */
        template<class T_SUBVEC> void set(const std::vector<PetscInt> &e_idxs,
                                          const T_SUBVEC &sub_vec)
        {
            VecSetValues(_v, e_idxs.size(), &e_idxs[0], &sub_vec[0], INSERT_VALUES);
        }

        /*!
           Get several entries
           \param e_idxs  Indicies of entries to be gotten.
                          Note: size_t cannot be the type of e_idxs template argument
           \param sub_vec Values of entries
        */
        template<class T_SUBVEC> void get(const std::vector<PetscInt> &e_idxs,
                                          T_SUBVEC &sub_vec)
        {
            VecGetValues(_v, e_idxs.size(), &e_idxs[0], &sub_vec[0]);
        }

        /*!
           Get local vector, i.e. entries in the same rank

           \param loc_vec  Pointer to array where stores the local vector,
                           memory allocation is not needed
        */
        PetscScalar *getLocalVector() const
        {
            PetscScalar *loc_vec;
            VecGetArray(_v, &loc_vec);
            return loc_vec;
        }

        /*!
           Get global vector

           \param u Array to store the global vector. Memory allocation is needed in advance
        */
        void getGlobalVector(PetscScalar u[]);

        /// Get an entry value. This is an expensive operation,
        /// and it only get local value. Use it for only test purpose
        PetscScalar get(const PetscInt idx) const
        {
            double x;
            VecGetValues(_v, 1, &idx, &x);
            return x;
        }

        /// Get PETsc vector. Use it only for test purpose
        PETSc_Vec &getData()
        {
            return _v;
        }

        /// Initialize the vector with a constant value
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
            \param file_name  File name for output
            \vw_format        File format listed as:
             PETSC_VIEWER_DEFAULT            Default format
             PETSC_VIEWER_ASCII_MATLAB       MATLAB format
             PETSC_VIEWER_ASCII_DENSE        Print matrix as dense
             PETSC_VIEWER_ASCII_IMPL         Implementation-specific format
                                               (which is in many cases the same as the default)
             PETSC_VIEWER_ASCII_INFO         Basic information about object
             PETSC_VIEWER_ASCII_INFO_DETAIL  More detailed info about object
             PETSC_VIEWER_ASCII_COMMON       Identical output format for all objects of a particular type
             PETSC_VIEWER_ASCII_INDEX        (for vectors) Prints the vector element number next to
                                                each vector entry
             PETSC_VIEWER_ASCII_SYMMODU      Print parallel vectors without indicating the processor ranges
             PETSC_VIEWER_ASCII_VTK          Outputs the object to a VTK file
             PETSC_VIEWER_NATIVE             Store the object to the binary file in its native format
                                              (for example, dense matrices are stored as dense),
                                              DMDA vectors are dumped directly to the file instead of
                                              being first put in the natural ordering
             PETSC_VIEWER_DRAW_BASIC         Views the vector with a simple 1d plot
             PETSC_VIEWER_DRAW_LG            Views the vector with a line graph
             PETSC_VIEWER_DRAW_CONTOUR       Views the vector with a contour plot
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

        /*!
              \brief  Collect local vectors
              \param  local_array Local array
              \param  global_array Global array
        */
        void gatherLocalVectors(PetscScalar local_array[],
                                PetscScalar global_array[]);

        friend void finalizeVectorAssembly(PETScVector &vec);
};

void finalizeVectorAssembly(PETScVector &vec);

} // end namespace
#endif
