/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <filesystem>

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace GeoLib
{
struct NamedRaster;
struct MinMaxPoints;
namespace IO
{
GeoLib::NamedRaster readRaster(BaseLib::ConfigTree const& raster_config,
                               std::string const& raster_directory,
                               GeoLib::MinMaxPoints const& min_max_points);
}  // namespace IO
}  // namespace GeoLib
