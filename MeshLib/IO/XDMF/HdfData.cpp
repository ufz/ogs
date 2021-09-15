/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "HdfData.h"

#include <hdf5.h>

#include <map>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "partition.h"

namespace MeshLib::IO
{
static hid_t meshPropertyType2HdfType(MeshPropertyDataType const ogs_data_type)
{
    std::map<MeshPropertyDataType const, hid_t> ogs_to_hdf_type = {
        {MeshPropertyDataType::float64, H5T_NATIVE_DOUBLE},
        {MeshPropertyDataType::float32, H5T_NATIVE_FLOAT},
        {MeshPropertyDataType::int32, H5T_NATIVE_INT32},
        {MeshPropertyDataType::int64, H5T_NATIVE_INT64},
        {MeshPropertyDataType::uint32, H5T_NATIVE_UINT32},
        {MeshPropertyDataType::uint64, H5T_NATIVE_UINT64},
        {MeshPropertyDataType::int8, H5T_NATIVE_INT8},
        {MeshPropertyDataType::uint8, H5T_NATIVE_UINT8}};
    try
    {
        return ogs_to_hdf_type.at(ogs_data_type);
    }
    catch (std::exception const& e)
    {
        OGS_FATAL("No known HDF5 type for OGS type. {:s}", e.what());
    }
}

HdfData::HdfData(void const* data_start, std::size_t const size_partitioned_dim,
                 std::size_t const size_tuple, std::string const& name,
                 MeshPropertyDataType const mesh_property_data_type,
                 unsigned int const n_files)
    : data_start(data_start), name(name)
{
    auto const& partition_info =
        getPartitionInfo(size_partitioned_dim, n_files);
    auto const& offset_partitioned_dim = partition_info.local_offset;
    offsets = {offset_partitioned_dim, 0};

    std::size_t unified_length = partition_info.local_length;

    chunk_space =
        (size_tuple > 1)
            ? std::vector<Hdf5DimType>{partition_info.longest_local_length,
                                       size_tuple}
            : std::vector<Hdf5DimType>{partition_info.longest_local_length};

    data_space = (size_tuple > 1)
                     ? std::vector<Hdf5DimType>{unified_length, size_tuple}
                     : std::vector<Hdf5DimType>{unified_length};
    file_space =
        (size_tuple > 1)
            ? std::vector<Hdf5DimType>{partition_info.global_length, size_tuple}
            : std::vector<Hdf5DimType>{partition_info.global_length};

    data_type = meshPropertyType2HdfType(mesh_property_data_type);

    DBUG(
        "HDF: dataset name: {:s}, offset: {:d}, data_space: {:d}, chunk_space "
        "{:d}, file_space: {:d}, tuples: {:d}",
        name, partition_info.local_offset, data_space[0], chunk_space[0],
        file_space[0], size_tuple);
}
}  // namespace MeshLib::IO