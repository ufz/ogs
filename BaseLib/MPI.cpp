// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MPI.h"

#ifdef USE_PETSC
namespace BaseLib::MPI
{
MPI_Comm OGS_COMM_WORLD = MPI_COMM_WORLD;
}
#endif  // USE_PETSC
