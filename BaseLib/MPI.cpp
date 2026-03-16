// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MPI.h"

#ifdef USE_PETSC
namespace BaseLib::MPI
{
MPI_Comm OGS_COMM_WORLD = MPI_COMM_WORLD;
}
#endif  // USE_PETSC

namespace BaseLib::MPI
{
Setup::Setup(int argc, char* argv[])
{
#ifdef USE_PETSC
    {
        int mpi_init;
        MPI_Initialized(&mpi_init);
        if (mpi_init == 1)
        {
            OGS_FATAL(
                "MPI has already been initialized. OGS does not support "
                "multiple MPI sessions in the same process");
        }
    }

    {
        int mpi_fin;
        MPI_Finalized(&mpi_fin);
        if (mpi_fin == 1)
        {
            OGS_FATAL(
                "MPI has already been shut down. OGS does not support multiple "
                "MPI sessions in the same process");
        }
    }

    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif
}

Setup::~Setup()
{
#ifdef USE_PETSC
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 1)
    {
        MPI_Finalize();
    }
#endif  // USE_PETSC
}
}  // namespace BaseLib::MPI
