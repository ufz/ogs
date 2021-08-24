/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  Collects and holds all metadata for writing XDMF file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <boost/shared_ptr.hpp>
#include <optional>
#include <string>
#include <vector>

#include "MeshLib/Location.h"
#include "MeshPropertyDataType.h"

class XdmfArrayType;

namespace MeshLib::IO
{
using XdmfDimType = unsigned int;

struct XdmfData final
{
    /**
     * \brief XdmfData contains meta data to be passed to the XdmfWriter - it
     * does not contain the actual values!!
     * @param size_partitioned_dim The first dimension (index 0) is assumed to
     * be the dimension that is partitioned. This first dimension expresses  the
     * length, usually the number of nodes or the number of cells. These values
     * give the length of the local partition.
     * @param size_tuple We assume there is at most a rank of 2 of data
     * (properties). The size of tuple gives the length of the second dimension
     * (index 1).
     * @param mesh_property_data_type property vector data type.
     * @param name The name of the attribute. It assumed to be unique.
     * @param attribute_center XdmfData is used for topology, geometry and
     * attributes. Geometry and topology have never a attribute_center.
     * Attributes have always an  attribute_center
     * @param index The position of the DataItem parents in a grid
     * (representing a single step). Convention is: 1=Time, 2=
     * Geometry, 3=Topology, 4>=Attribute
     * @param num_of_files If greater than 1 it groups the data of each process.
     * The num_of_files specifies the number of groups
     *
     */
    XdmfData(std::size_t size_partitioned_dim, std::size_t size_tuple,
             MeshPropertyDataType mesh_property_data_type,
             std::string const& name,
             std::optional<MeshLib::MeshItemType> attribute_center,
             unsigned int const index, unsigned int num_of_files);
    // a hyperslab is defined by starts and strides see
    // https://www.xdmf.org/index.php/XDMF_Model_and_Format#HyperSlab
    std::vector<XdmfDimType> starts;
    std::vector<XdmfDimType> strides;
    std::vector<XdmfDimType> global_block_dims;
    MeshPropertyDataType data_type;
    std::string name;
    std::optional<MeshLib::MeshItemType> attribute_center;
    unsigned int index;
};

}  // namespace MeshLib::IO
