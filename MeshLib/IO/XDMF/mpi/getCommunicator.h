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
FileCommunicator getCommunicator(unsigned int num_of_files);
}  // namespace MeshLib::IO