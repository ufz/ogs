/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <algorithm>

#include "Error.h"

#ifdef USE_PETSC
#include <mpi.h>
#endif

namespace BaseLib::MPI
{

#ifdef USE_PETSC
struct Mpi
{
    Mpi(MPI_Comm const communicator = MPI_COMM_WORLD)
        : communicator(communicator)
    {
        int mpi_init;
        MPI_Initialized(&mpi_init);
        if (mpi_init != 1)
        {
            OGS_FATAL("MPI is not initialized.");
        }
        MPI_Comm_size(communicator, &size);
        MPI_Comm_rank(communicator, &rank);
    }

    MPI_Comm communicator;
    int size;
    int rank;
};

template <typename T>
constexpr MPI_Datatype mpiType()
{
    using U = std::remove_const_t<T>;
    if constexpr (std::is_same_v<U, bool>)
    {
        return MPI_C_BOOL;
    }
    if constexpr (std::is_same_v<U, char>)
    {
        return MPI_CHAR;
    }
    if constexpr (std::is_same_v<U, double>)
    {
        return MPI_DOUBLE;
    }
    if constexpr (std::is_same_v<U, float>)
    {
        return MPI_FLOAT;
    }
    if constexpr (std::is_same_v<U, int>)
    {
        return MPI_INT;
    }
    if constexpr (std::is_same_v<U, std::size_t>)
    {
        return MPI_UNSIGNED_LONG;
    }
    if constexpr (std::is_same_v<U, unsigned int>)
    {
        return MPI_UNSIGNED;
    }
}

template <typename T>
static std::vector<T> allgather(T const& value, Mpi const& mpi)
{
    std::vector<T> result(mpi.size);

    result[mpi.rank] = value;

    MPI_Allgather(&result[mpi.rank], 1, mpiType<T>(), result.data(), 1,
                  mpiType<T>(), mpi.communicator);

    return result;
}

template <typename T>
static std::vector<T> allgather(std::vector<T> const& vector, Mpi const& mpi)
{
    std::size_t const size = vector.size();
    // Flat in memory over all ranks;
    std::vector<T> result(mpi.size * size);

    std::copy_n(vector.begin(), size, &result[mpi.rank * size]);

    MPI_Allgather(&result[mpi.rank * size], size, mpiType<T>(), result.data(),
                  size, mpiType<T>(), mpi.communicator);

    return result;
}

template <typename T>
static T allreduce(T const& value, MPI_Op const& mpi_op, Mpi const& mpi)
{
    T result{};

    MPI_Allreduce(&value, &result, 1, mpiType<T>(), mpi_op, mpi.communicator);
    return result;
}

template <typename T>
static void allreduceInplace(std::vector<T>& vector,
                             MPI_Op const& mpi_op,
                             Mpi const& mpi)
{
    MPI_Allreduce(MPI_IN_PLACE,
                  vector.data(),
                  vector.size(),
                  mpiType<T>(),
                  mpi_op,
                  mpi.communicator);
}
#endif

/// The reduction is implemented transparently for with and without MPI. In the
/// latter case the input value is returned.
static inline int reduceMin(int const val)
{
#ifdef USE_PETSC
    return allreduce(val, MPI_MIN, Mpi{MPI_COMM_WORLD});
#else
    return val;
#endif
}

}  // namespace BaseLib::MPI
