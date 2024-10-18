/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Error.h"

#ifdef USE_PETSC
#include <mpi.h>
#endif

namespace BaseLib::MPI
{

static inline int reduceMin(int const val)
{
#ifdef USE_PETSC
    // Reduce operations for interprocess communications while using Petsc
    int result;
    MPI_Allreduce(&val, &result, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD);
    return result;
#else
    // Reduce operations for interprocess communications without using Petsc
    return val;
#endif
}

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
#endif
}  // namespace BaseLib::MPI
