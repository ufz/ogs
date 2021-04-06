/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

std::pair<std::size_t, std::size_t> getPartitionInfo(std::size_t const size)
{
    return {0, size};
}
}  // namespace MeshLib::IO
