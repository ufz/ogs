// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <mpi.h>

#include <filesystem>

namespace MeshLib::IO
{
struct FileCommunicator final
{
    MPI_Comm mpi_communicator;
    int color;
    std::filesystem::path output_filename;
};
FileCommunicator getCommunicator(unsigned int n_files);
}  // namespace MeshLib::IO