/**
 * \file
 * \author Tobias Meisel
 * \date   2021-07-13
 * \brief  write_xdmf generates a function based on spatial mesh data. The
 * generated function finally generates an XDMF string when temporal data is
 * known
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "MeshLib/IO/XDMF/XdmfData.h"
#include "MeshLib/Location.h"
#include "MeshPropertyDataType.h"

using XdmfDimType = unsigned int;

namespace MeshLib::IO
{
/**
 * \brief Generator function that creates a function capturing the spatial data
 * of a mesh Temporal data can later be passed as argument
 * @param geometry Metadata for the geometry (points) of the mesh
 * @param topology Metadata for the topology of the mesh
 * @param variable_attributes Meta data for attributes changing over time
 * @param constant_attributes Meta data for attributes NOT changing over time
 * @param h5filename Name of the file where the actual data was written
 * @param ogs_version OGS Version to be added to XdmfInformation tag
 * @param mesh_name Name of the output mesh
 * @return unary function with vector of time step values, returning XDMF string
 */
std::function<std::string(std::vector<double>)> write_xdmf(
    XdmfData const& geometry, XdmfData const& topology,
    std::vector<XdmfData> const& constant_attributes,
    std::vector<XdmfData> const& variable_attributes,
    std::string const& h5filename, std::string const& ogs_version,
    std::string const& mesh_name);
}  // namespace MeshLib::IO