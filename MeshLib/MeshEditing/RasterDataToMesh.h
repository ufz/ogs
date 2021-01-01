/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GeoLib/Raster.h"
#include "MeshLib/Mesh.h"


namespace MeshLib
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
}  // end namespace MeshLib
