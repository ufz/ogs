/**
 * \file
 * \brief Definition of member functions of class PETScVector,
 *        which provides an interface to PETSc vector routines.
 *
 *   Note: the return message of PETSc routines is omitted in
 *         the source code. If it is really needed, it can be activated by
 *         adding a PetscErrorCode type variable before each PETSc function
 *
 * \author Wenqing Wang
 * \date Nov 2011 - Sep 2013
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PETScVector.h"

#include <algorithm>
#include <cassert>

#include "BaseLib/Error.h"

namespace MathLib
{
PETScVector::PETScVector(const PetscInt vec_size, const bool is_global_size)
{
    if (is_global_size)
    {
        VecCreate(PETSC_COMM_WORLD, &v_);
        VecSetSizes(v_, PETSC_DECIDE, vec_size);
    }
    else
    {
        // Fix size partitioning
        // the size can be associated to specific memory allocation of a matrix
        VecCreateMPI(PETSC_COMM_WORLD, vec_size, PETSC_DECIDE, &v_);
    }

    config();
}

PETScVector::PETScVector(const PetscInt vec_size,
                         const std::vector<PetscInt>& ghost_ids,
                         const bool is_global_size)
    : size_ghosts_{static_cast<PetscInt>(ghost_ids.size())},
      created_with_ghost_id_{true}
{
    PetscInt nghosts = static_cast<PetscInt>(ghost_ids.size());
    if (is_global_size)
    {
        VecCreateGhost(PETSC_COMM_WORLD, PETSC_DECIDE, vec_size, nghosts,
                       ghost_ids.data(), &v_);
    }
    else
    {
        VecCreate(PETSC_COMM_WORLD, &v_);
        VecSetType(v_, VECMPI);
        VecSetSizes(v_, vec_size, PETSC_DECIDE);
        VecMPISetGhost(v_, nghosts, ghost_ids.data());
    }

    config();

    for (PetscInt i = 0; i < nghosts; i++)
    {
        global_ids2local_ids_ghost_.emplace(ghost_ids[i], size_loc_ + i);
    }
}

PETScVector::PETScVector(const PETScVector& existing_vec, const bool deep_copy)
{
    shallowCopy(existing_vec);

    // Copy values
    if (deep_copy)
    {
        VecCopy(existing_vec.v_, v_);
    }
}

PETScVector::PETScVector(PETScVector&& other)
    : v_{std::move(other.v_)},
      v_loc_{std::move(other.v_loc_)},
      start_rank_{other.start_rank_},
      end_rank_{other.end_rank_},
      size_{other.size_},
      size_loc_{other.size_loc_},
      size_ghosts_{other.size_ghosts_},
      created_with_ghost_id_{other.created_with_ghost_id_},
      global_ids2local_ids_ghost_{other.global_ids2local_ids_ghost_}
{
}

void PETScVector::config()
{
    VecSetFromOptions(v_);
    VecGetOwnershipRange(v_, &start_rank_, &end_rank_);

    VecGetLocalSize(v_, &size_loc_);
    VecGetSize(v_, &size_);

    VecSetOption(v_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}

void PETScVector::finalizeAssembly()
{
    VecAssemblyBegin(v_);
    VecAssemblyEnd(v_);
}

void PETScVector::gatherLocalVectors(PetscScalar local_array[],
                                     PetscScalar global_array[]) const
{
    // Collect vectors from processors.
    int size_rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size_rank);

    // number of elements to be sent for each rank
    std::vector<PetscInt> i_cnt(size_rank);

    MPI_Allgather(&size_loc_, 1, MPI_INT, &i_cnt[0], 1, MPI_INT,
                  PETSC_COMM_WORLD);

    // collect local array
    PetscInt offset = 0;
    // offset in the receive vector of the data from each rank
    std::vector<PetscInt> i_disp(size_rank);
    for (PetscInt i = 0; i < size_rank; i++)
    {
        i_disp[i] = offset;
        offset += i_cnt[i];
    }

    MPI_Allgatherv(local_array, size_loc_, MPI_DOUBLE, global_array, &i_cnt[0],
                   &i_disp[0], MPI_DOUBLE, PETSC_COMM_WORLD);
}

void PETScVector::getGlobalVector(std::vector<PetscScalar>& u) const
{
#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    assert(static_cast<PetscInt>(u.size()) == size_);

    PetscInt state;
    VecLockGet(v_, &state);
    if (state != 0)
    {
        OGS_FATAL("PETSc vector is already locked for {:s} access.",
                  state > 0 ? "read" : "write");
    }
    PetscScalar* xp = nullptr;
    VecGetArray(v_, &xp);

    gatherLocalVectors(xp, u.data());

    // This following line may be needed late on
    //  for a communication load balance:
    // MPI_Barrier(PETSC_COMM_WORLD);

    VecRestoreArray(v_, &xp);

// TEST
#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(
        PETSC_COMM_WORLD,
        "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1,
        mem2, (int)(mem2 - mem1));
#endif
}

void PETScVector::setLocalAccessibleVector() const
{
    if (created_with_ghost_id_)
    {
        copyValues(entry_array_);
        return;
    }

    entry_array_.resize(size_);
    getGlobalVector(entry_array_);
}

void PETScVector::copyValues(std::vector<PetscScalar>& u) const
{
    u.resize(getLocalSize() + getGhostSize());

    PetscScalar* loc_x = getLocalVector();
    std::copy_n(loc_x, getLocalSize() + getGhostSize(), u.begin());
    restoreArray(loc_x);
}

PetscScalar PETScVector::get(const PetscInt idx) const
{
    if (created_with_ghost_id_)
    {
        return entry_array_[getLocalIndex(idx)];
    }

    return entry_array_[idx];
}

std::vector<PetscScalar> PETScVector::get(
    std::vector<IndexType> const& indices) const
{
    std::vector<PetscScalar> local_x(indices.size());
    std::transform(indices.begin(), indices.end(), local_x.begin(),
                   [this](IndexType index) { return get(index); });
    return local_x;
}

PetscScalar* PETScVector::getLocalVector() const
{
    PetscScalar* loc_array;
    if (created_with_ghost_id_ && !global_ids2local_ids_ghost_.empty())
    {
        VecGhostUpdateBegin(v_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(v_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostGetLocalForm(v_, &v_loc_);
        VecGetArray(v_loc_, &loc_array);
    }
    else
    {
        VecGetArray(v_, &loc_array);
    }

    return loc_array;
}

PetscInt PETScVector::getLocalIndex(const PetscInt global_index) const
{
    if (global_index >= 0)  // non-ghost entry.
    {
#ifndef NDEBUG
        if (global_index < start_rank_ || global_index >= end_rank_)
        {
            OGS_FATAL(
                "The global index {:d} is out of the range `[`{:d}, {:d}`)` of "
                "the current rank.",
                global_index, start_rank_, end_rank_);
        }
#endif
        return global_index - start_rank_;
    }

    // A special case for a ghost location with global index equal to
    // the size of the local vector:
    PetscInt real_global_index = (-global_index == size_) ? 0 : -global_index;

#ifndef NDEBUG
    if (global_ids2local_ids_ghost_.find(real_global_index) ==
            global_ids2local_ids_ghost_.end() ||
        global_ids2local_ids_ghost_.empty())
    {
        OGS_FATAL("The global index {:d} is not found as a ghost ID",
                  global_index);
    }
#endif

    return global_ids2local_ids_ghost_.at(real_global_index);
}

void PETScVector::restoreArray(PetscScalar* array) const
{
    if (created_with_ghost_id_ && !global_ids2local_ids_ghost_.empty())
    {
        VecRestoreArray(v_loc_, &array);
        VecGhostRestoreLocalForm(v_, &v_loc_);
    }
    else
    {
        VecRestoreArray(v_, &array);
    }
}

void PETScVector::viewer(const std::string& file_name,
                         const PetscViewerFormat vw_format) const
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    PetscObjectSetName((PetscObject)v_, file_name.c_str());
    VecView(v_, viewer);

#define nEXIT_TEST
#ifdef EXIT_TEST
    VecDestroy(v_);
    PetscFinalize();
    exit(0);
#endif
}

void PETScVector::shallowCopy(const PETScVector& v)
{
    destroy();

    VecDuplicate(v.getRawVector(), &v_);

    start_rank_ = v.start_rank_;
    end_rank_ = v.end_rank_;
    size_ = v.size_;
    size_loc_ = v.size_loc_;
    size_ghosts_ = v.size_ghosts_;
    created_with_ghost_id_ = v.created_with_ghost_id_;
    global_ids2local_ids_ghost_ = v.global_ids2local_ids_ghost_;

    config();
}

void finalizeVectorAssembly(PETScVector& vec)
{
    vec.finalizeAssembly();
}

}  // namespace MathLib
