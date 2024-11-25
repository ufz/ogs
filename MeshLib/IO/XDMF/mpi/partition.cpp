/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Function specific to execution with MPI!!
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "../partition.h"

#include <mpi.h>

#include <numeric>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "MeshLib/IO/XDMF/fileIO.h"
#include "getCommunicator.h"

namespace MeshLib::IO
{
bool isFileManager()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    return mpi_rank == 0;
}

PartitionInfo getPartitionInfo(std::size_t const size,
                               unsigned int const n_files)
{
    BaseLib::MPI::Mpi const mpi{getCommunicator(n_files).mpi_communicator};

    std::vector<std::size_t> const partition_sizes =
        BaseLib::MPI::allgather(size, mpi);

    // the first partition's offset is zero, offsets for subsequent
    // partitions are the accumulated sum of all preceding size.
    std::vector<std::size_t> const partition_offsets =
        BaseLib::sizesToOffsets(partition_sizes);

    // chunked
    std::size_t longest_partition =
        *max_element(partition_sizes.begin(), partition_sizes.end());

    // local_offset, local_length, longest_local_length, global_length
    return {partition_offsets[mpi.rank], size, longest_partition,
            partition_offsets.back()};
}
}  // namespace MeshLib::IO
