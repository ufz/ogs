// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "MeshLib/IO/XDMF/XdmfData.h"
#include "MeshPropertyDataType.h"

using XdmfDimType = unsigned int;

namespace MeshLib::IO
{
/**
 * \brief Generator function that creates a function capturing the spatial data
 * of a mesh Temporal data can later be passed as argument
 * \param geometry Metadata for the geometry (points) of the mesh
 * \param topology Metadata for the topology of the mesh
 * \param variable_attributes Meta data for attributes changing over time
 * \param constant_attributes Meta data for attributes NOT changing over time
 * \param h5filename Name of the file where the actual data was written
 * \param ogs_version OGS Version to be added to XdmfInformation tag
 * \param mesh_name Name of the output mesh
 * \return unary function with vector of time step values, returning XDMF string
 */
std::function<std::string(std::vector<double>)> write_xdmf(
    XdmfData const& geometry, XdmfData const& topology,
    std::vector<XdmfData> const& constant_attributes,
    std::vector<XdmfData> const& variable_attributes,
    std::string const& h5filename, std::string const& ogs_version,
    std::string const& mesh_name);
}  // namespace MeshLib::IO
