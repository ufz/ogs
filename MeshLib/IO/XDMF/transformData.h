/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  Transforms OGS Mesh into vectorized data.
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "XdmfHdfData.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshLib::IO
{
/**
 * \brief Create meta data for attributes used for hdf5 and xdmf
 * @param mesh OGS mesh can be mesh or partitionedMesh
 * @return vector of meta data
 */
std::vector<AttributeMeta> transformAttributes(MeshLib::Mesh const& mesh);
/**
 * \brief Copies all node points into a new vector. Contiguous data used for
 * writing.
 * @param mesh OGS mesh can be mesh or partitionedMesh
 * @return Geometry containing a copy of the data and the geometry meta data
 */
Geometry transformGeometry(MeshLib::Mesh const& mesh);
/**
 * \brief Copies all cells into a new vector. Contiguous data used for writing.
 * The topology is specific to xdmf because it contains the xdmf cell types!!
 * See section topology in https://www.xdmf.org/index.php/XDMF_Model_and_Format
 * @param mesh OGS mesh can be mesh or partitionedMesh
 * @return Topology containing a copy of the data and the topology meta data
 */
Topology transformTopology(MeshLib::Mesh const& mesh, std::size_t offset);
}  // namespace MeshLib::IO