/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Collects and holds all metadata for writing HDF5 file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>
#include <vector>

#include "MeshPropertyDataType.h"

namespace MeshLib::IO
{
using Hdf5DimType = unsigned long long;

struct HdfData final
{
    HdfData(void const* data_start, std::size_t size_partitioned_dim,
            std::size_t size_tuple, std::string const& name,
            MeshPropertyDataType mesh_property_data_type);
    void const* data_start;
    std::vector<Hdf5DimType> data_space;
    std::vector<Hdf5DimType> offsets;
    std::vector<Hdf5DimType> file_space;
    std::vector<Hdf5DimType> chunk_space;
    std::string name;
    int64_t data_type;
};

}  // namespace MeshLib::IO