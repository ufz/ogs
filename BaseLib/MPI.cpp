/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MPI.h"

#ifdef USE_PETSC
namespace BaseLib::MPI
{
MPI_Comm OGS_COMM_WORLD = MPI_COMM_WORLD;
}
#endif  // USE_PETSC
