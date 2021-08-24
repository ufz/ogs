/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Function specific to execution with MPI, never include directly!!
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "getCommunicator.h"

#include <hdf5.h>
#include <math.h>
#include <mpi.h>

#include <cassert>

#include "BaseLib/Logging.h"
#include "MeshLib/IO/XDMF/fileIO.h"

using namespace std::string_literals;

namespace MeshLib::IO
{
int getGroupIndex(int const input_index, int const input_size,
                  int const new_group_size)
{
    // A grouping algorithm that determines the number of groups and return the
    // group idx of the specified input_index
    assert(input_size >= new_group_size);
    int const minimal_output_group_size =
        std::lround(input_size / new_group_size);
    int const maximal_output_group_size = (input_size % new_group_size)
                                              ? minimal_output_group_size + 1
                                              : minimal_output_group_size;
    return std::lround(input_index / maximal_output_group_size);
};

FileCommunicator getCommunicator(unsigned int const num_of_files)
{
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    int rank_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
    int file_group_id = getGroupIndex(rank_id, num_procs, num_of_files);
    MPI_Comm new_communicator;
    MPI_Comm_split(MPI_COMM_WORLD, file_group_id, rank_id, &new_communicator);
    return FileCommunicator{std::move(new_communicator),
                            std::move(file_group_id), ""};
}
}  // namespace MeshLib::IO