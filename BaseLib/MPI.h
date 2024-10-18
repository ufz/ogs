/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

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
    MPI_Allreduce(&val, &result, 1, MPI_INTEGER, MPI_MIN, PETSC_COMM_WORLD);
    return result;
#else
    // Reduce operations for interprocess communications without using Petsc
    return val;
#endif
}

}  // namespace BaseLib::MPI
