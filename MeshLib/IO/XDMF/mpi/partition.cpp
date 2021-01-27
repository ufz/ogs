/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Function specific to execution with MPI!!
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <mpi.h>

#include <deque>
#include <numeric>


#include "../partition.h"
#include "BaseLib/Logging.h"
#include "MeshLib/IO/XDMF/fileIO.h"

namespace MeshLib::IO
{
bool isFileManager()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    return mpi_rank == 0;
}

PartitionInfo getPartitionInfo(std::size_t const size)
{
    MPI_Comm const mpi_comm = MPI_COMM_WORLD;
    int mpi_size;
    int mpi_rank;
    MPI_Comm_size(mpi_comm, &mpi_size);
    MPI_Comm_rank(mpi_comm, &mpi_rank);

    std::vector<std::size_t> partition_sizes;
    partition_sizes.resize(mpi_size);

    MPI_Allgather(&size,
                  1,
                  MPI_UNSIGNED_LONG,
                  partition_sizes.data(),
                  1,
                  MPI_UNSIGNED_LONG,
                  mpi_comm);

    // the first partition's offset is zero, offsets for subsequent
    // partitions are the accumulated sum of all preceding size (excluding
    // own size)
    std::vector<std::size_t> partition_offsets(1, 0);
    std::partial_sum(partition_sizes.begin(),
                     partition_sizes.end(),
                     back_inserter(partition_offsets));

    //chunked
    std::size_t longest_partition = *max_element(partition_sizes.begin(), partition_sizes.end());
    auto this_partition_offset = longest_partition*mpi_rank;
    return {this_partition_offset, longest_partition,partition_sizes.size()};
    //return {partition_offsets[mpi_rank], partition_offsets.back()};
}
}  // namespace MeshLib::IO
