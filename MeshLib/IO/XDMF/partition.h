/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Dispatches functions specific to execution platform (w/o MPI)
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <utility>

namespace MeshLib::IO
{
struct PartitionInfo
{
    std::size_t local_offset;
    std::size_t local_length;
    std::size_t longest_local_length;
    std::size_t global_length;
};

PartitionInfo getPartitionInfo(std::size_t const size,
                               unsigned int const n_files);
bool isFileManager();
}  // namespace MeshLib::IO
