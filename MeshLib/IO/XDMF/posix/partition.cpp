/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshLib/IO/XDMF/partition.h"

namespace MeshLib::IO
{
bool isFileManager()
{
    return true;
}

PartitionInfo getPartitionInfo(std::size_t const size, unsigned int)
{
    // local_offset, local_length, longest_local_length, global_number_process
    return {0, size, size, size};
}
}  // namespace MeshLib::IO
