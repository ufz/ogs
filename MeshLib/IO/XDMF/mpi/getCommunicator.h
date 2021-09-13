/**
 * \file
 * \author Tobias Meisel
 * \date   2021-09-14
 * \brief  Assigns to each MPI communicator an output file name by attribute
 * color There are multiple implementation to this interface!
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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