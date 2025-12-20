// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <hdf5.h>

#include <string>
#include <vector>

#include "MeshPropertyDataType.h"

namespace MeshLib::IO
{
using Hdf5DimType = hsize_t;

struct HdfData final
{
    HdfData(void const* data_start, std::size_t size_partitioned_dim,
            std::size_t size_tuple, std::string const& name,
            MeshPropertyDataType mesh_property_data_type, unsigned int n_files,
            unsigned int chunk_size_bytes);
    void const* data_start;
    std::vector<Hdf5DimType> data_space;
    std::vector<Hdf5DimType> offsets;
    std::vector<Hdf5DimType> file_space;
    std::vector<Hdf5DimType> chunk_space;
    std::string name;
    int64_t data_type;
};

}  // namespace MeshLib::IO
