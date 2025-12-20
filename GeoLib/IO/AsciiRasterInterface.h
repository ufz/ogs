// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "GeoLib/Raster.h"

namespace FileIO
{
/**
 * Interface for reading and writing a number of ASCII raster formats.
 * Currently supported are reading and writing of Esri asc-files and
 * reading of Surfer grd-files.
 */
class AsciiRasterInterface
{
public:
    /// Reads raster file by detecting type based on extension and then calling
    /// the appropriate method
    static GeoLib::Raster* readRaster(std::string const& fname);

    /// Reads an ArcGis ASC raster file
    static GeoLib::Raster* getRasterFromASCFile(std::string const& fname);

    /// Reads a Surfer GRD raster file
    static GeoLib::Raster* getRasterFromSurferFile(std::string const& fname);

    /// Reads a XYZ raster file
    static GeoLib::Raster* getRasterFromXyzFile(std::string const& fname);

    /// Writes an Esri asc-file
    static void writeRasterAsASC(GeoLib::Raster const& raster,
                                 std::string const& file_name);
};

/// Reads a vector of rasters given by file names. On error nothing is returned,
/// otherwise the returned vector contains pointers to the read rasters.
std::optional<std::vector<GeoLib::Raster const*>> readRasters(
    std::vector<std::string> const& raster_paths);
}  // end namespace FileIO
