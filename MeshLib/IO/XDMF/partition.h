// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
