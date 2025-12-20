// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
