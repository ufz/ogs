/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "XdmfData.h"

#include <cassert>
#include <map>

#include "BaseLib/Error.h"
#include "MeshLib/Location.h"
#include "partition.h"

namespace MeshLib::IO
{
XdmfData::XdmfData(std::size_t const size_partitioned_dim,
                   std::size_t const size_tuple,
                   MeshPropertyDataType const mesh_property_data_type,
                   std::string const& name,
                   std::optional<MeshLib::MeshItemType> const attribute_center,
                   unsigned int const index)
    : starts([&size_tuple]() {
          if (size_tuple > 1)
          {
              return std::vector<XdmfDimType>{0, 0};
          }
          else
          {
              return std::vector<XdmfDimType>{0};
          }
      }()),
      strides([&size_tuple]() {
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
      name(name),
      attribute_center(attribute_center),
      index(index)
{
    auto partition_info = getPartitionInfo(size_partitioned_dim);
    // TODO (tm) XdmfLib does not support 64 bit data types so far
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