/**
 * \file  PETScVector.h
 * \brief Declaration of class PETScVector, which provides an interface to
 *        PETSc vector routines.
 *
 * \author Wenqing Wang
 * \date Nov 2011 - Sep 2013
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include <petscvec.h>

namespace MathLib
{
/*!
   \class PETScVector

   \brief Wrapper class for PETSc vector
*/
class PETScVector
{
    public:
        using IndexType = PetscInt;
        // TODO make this class opaque, s.t. the PETSc symbols are not scattered all
        //      over the global namespace
        using PETSc_Vec = Vec;

    public:
        PETScVector() {}

        /*!
            \brief Constructor
            \param vec_size       The size of the vector, either global or local
            \param is_global_size The flag of the type of vec_size, i.e.
                                  whether it is a global size
                                  or local size. The default is true.
                                  If is_global_size is true, the vector
                                  is created by the global size, the local size
                                  of the vector is determined by PETSc,
                                  and vice versa is the same.
        */
        PETScVector(const PetscInt vec_size, const bool is_global_size = true);

        /*!
            \brief Constructor
            \param vec_size       The size of the vector, either global or local
            \param ghost_ids      The global indices of ghost entries
            \param is_global_size The flag of the type of vec_size, i.e. whether it is a global size
                                  or local size. The default is true.
        */
        PETScVector(const PetscInt vec_size, const std::vector<PetscInt>& ghost_ids,
                    const bool is_global_size = true);

        /*!
             \brief Copy constructor
             \param existing_vec The vector to be copied
             \param deep_copy    The flag for a deep copy, which means to copy the values as well,
                                 the default is true
        */
        explicit PETScVector(const PETScVector &existing_vec, const bool deep_copy = true);

        PETScVector(PETScVector&& other);

        ~PETScVector()
        {
            destroy();
        }

        /// Perform MPI collection of assembled entries in buffer
        void finalizeAssembly();

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

        /// Get the number of ghost entries in the same rank
        PetscInt getGhostSize() const
        {
            return _size_ghosts;
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
                          Note: std::size_t cannot be the type of e_idxs template argument
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
                          Note: std::size_t cannot be the type of e_idxs template argument
           \param sub_vec Entries to be added
        */
        template<class T_SUBVEC> void set(const std::vector<PetscInt> &e_idxs,
                                          const T_SUBVEC &sub_vec)
        {
            VecSetValues(_v, e_idxs.size(), &e_idxs[0], &sub_vec[0], INSERT_VALUES);
        }

        // TODO preliminary
        void setZero() { VecSet(_v, 0.0); }

        /// Disallow moving.
        /// \todo This operator should be implemented properly when doing a
        ///       general cleanup of all matrix and vector classes.
        PETScVector& operator = (PETScVector &&) = delete;

        /// Set local accessible vector in order to get entries.
        /// Call this before call operator[] or get(...).
        void setLocalAccessibleVector() const;

        /// Get several entries. setLocalAccessibleVector() must be
        /// called beforehand.
        std::vector<double> get(std::vector<IndexType> const& indices) const;

        /// Get the value of an entry by [] operator.
        /// setLocalAccessibleVector() must be called beforehand.
        double operator[] (PetscInt idx) const
        {
            return get(idx);
        }

        /*!
           Get global vector
           \param u Array to store the global vector. Memory allocation is needed in advance
        */
        void getGlobalVector(std::vector<PetscScalar>& u) const;

        /* Get an entry value. This is an expensive operation,
           and it only get local value. Use it for only test purpose
           Get the value of an entry by [] operator.
           setLocalAccessibleVector() must be called beforehand.
        */
        PetscScalar get(const PetscInt idx) const;

        //! Exposes the underlying PETSc vector.
        PETSc_Vec getRawVector() { return _v; }

        /*! Exposes the underlying PETSc vector.
         *
         * \warning
         * This method is dangerous insofar as you can do arbitrary things also
         * with a const PETSc vector!
         */
        PETSc_Vec getRawVector() const { return _v; }

        /*!
           Copy local entries including ghost ones to an array
           \param u Preallocated vector for the values of local entries.
        */
        void copyValues(std::vector<double>& u) const;

        /*! View the global vector for test purpose. Do not use it for output a big vector.
            \param file_name  File name for output
            \param vw_format  File format listed as:
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
                    const PetscViewerFormat vw_format = PETSC_VIEWER_ASCII_MATLAB ) const;

        void shallowCopy(const PETScVector &v);

    private:
        void destroy() { if (_v) VecDestroy(&_v); _v = nullptr; }

        PETSc_Vec _v = nullptr;
        /// Local vector, which is only for the case that  _v is created
        /// with ghost entries.
        mutable PETSc_Vec _v_loc = nullptr;

        /// Starting index in a rank
        PetscInt _start_rank;
        /// Ending index in a rank
        PetscInt _end_rank;

        /// Size of the vector
        PetscInt _size;
        /// Size of local entries
        PetscInt _size_loc;
        /// Size of local ghost entries
        PetscInt _size_ghosts = 0;

        /// Flag to indicate whether the vector is created with ghost entry indices
        bool _has_ghost_id = false;

        /*!
           \brief Array containing the entries of the vector.
           If the vector is created without given ghost IDs, the array
           contains all entries of the global vector, _v. Otherwise it
           only contains the entries of the global vector owned by the
           current rank.
        */
        mutable std::vector<PetscScalar> _entry_array;

        /// Map global indices of ghost enrties to local indices
        mutable std::map<PetscInt, PetscInt> _global_ids2local_ids_ghost;

        /*!
              \brief  Collect local vectors
              \param  local_array Local array
              \param  global_array Global array
        */
        void gatherLocalVectors(PetscScalar local_array[],
                                PetscScalar global_array[]) const;

        /*!
           Get local vector, i.e. entries in the same rank
        */
        PetscScalar* getLocalVector() const;

        /// Get local index by a global index
        PetscInt getLocalIndex(const PetscInt global_index) const;

        /*!
           Restore array after finish access local array
           \param array  Pointer to the local array fetched by VecGetArray
        */
        inline void restoreArray(PetscScalar* array) const;

        /// A funtion called by constructors to configure members
        void config();

        friend void finalizeVectorAssembly(PETScVector &vec);
};

/// Function to finalize the vector assembly
void finalizeVectorAssembly(PETScVector &vec);

} // end namespace
