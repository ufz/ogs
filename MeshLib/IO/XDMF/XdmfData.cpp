// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "XdmfData.h"

#include <cassert>
#include <map>

#include "BaseLib/Error.h"
#include "partition.h"

namespace MeshLib::IO
{
XdmfData::XdmfData(std::size_t const size_partitioned_dim,
                   std::size_t const size_tuple,
                   MeshPropertyDataType const mesh_property_data_type,
                   std::string const& name,
                   std::optional<MeshLib::MeshItemType> const attribute_center,
                   unsigned int const index,
                   unsigned int const n_files,
                   std::optional<ParentDataType> const parent_data_type)
    : starts(
          [&size_tuple]()
          {
              if (size_tuple > 1)
              {
                  return std::vector<XdmfDimType>{0, 0};
              }
              else
              {
                  return std::vector<XdmfDimType>{0};
              }
          }()),
      strides(
          [&size_tuple]()
          {
              if (size_tuple > 1)
              {
                  return std::vector<XdmfDimType>{1, 1};
              }
              else
              {
                  return std::vector<XdmfDimType>{1};
              }
          }()),
      data_type(mesh_property_data_type),
      size_partitioned_dim(size_partitioned_dim),
      name(name),
      attribute_center(attribute_center),
      index(index),
      parent_data_type(parent_data_type)
{
    auto partition_info = getPartitionInfo(size_partitioned_dim, n_files);
    assert(partition_info.local_length <
           std::numeric_limits<unsigned int>::max());
    auto const ui_global_components =
        static_cast<unsigned int>(partition_info.global_length);
    auto const ui_tuple_size = static_cast<unsigned int>(size_tuple);

    if (ui_tuple_size == 1)
    {
        global_block_dims = {ui_global_components};
    }
    else
    {
        global_block_dims = {ui_global_components, ui_tuple_size};
    }

    DBUG(
        "XDMF: dataset name: {:s}, offset: {:d} "
        "global_blocks: {:d}, tuples: {:d}",
        name, partition_info.local_offset, global_block_dims[0], ui_tuple_size);
}
}  // namespace MeshLib::IO