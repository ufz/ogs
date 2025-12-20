// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "GeoLib/Raster.h"
#include "MeshLib/Mesh.h"

namespace MeshToolsLib
{

/**
 * \brief Adding pixel values from a raster onto nodes or cells of a mesh
 */
namespace RasterDataToMesh
{
bool projectToNodes(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                    double const default_replacement,
                    std::string const& array_name);

bool projectToElements(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                       double const default_replacement,
                       std::string const& array_name);
}  // end namespace RasterDataToMesh
}  // namespace MeshToolsLib
