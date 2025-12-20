// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
