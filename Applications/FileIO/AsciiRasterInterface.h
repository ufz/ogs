/**
 * @file AsciiRasterInterface.h
 * @author Karsten Rink
 * @date 2014-09-10
 * @brief Definition of the AsciiRasterInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef ASCIIRASTERINTERFACE_H_
#define ASCIIRASTERINTERFACE_H_

#include <fstream>
#include <vector>
#include <string>
#include <boost/optional.hpp>

#include "GeoLib/Raster.h"

namespace FileIO
{
/**
 * Interface for reading and writing a number of ASCII raster formats.
 * Currently supported are reading and writing of Esri asc-files and
 * reading of Surfer grd-files.
 */
class AsciiRasterInterface {
public:
    /// Reads raster file by detecting type based on extension and then calling the apropriate method
    static GeoLib::Raster* readRaster(std::string const& fname);

    /// Reads an ArcGis ASC raster file
    static GeoLib::Raster* getRasterFromASCFile(std::string const& fname);

    /// Reads a Surfer GRD raster file
    static GeoLib::Raster* getRasterFromSurferFile(std::string const& fname);

    /// Writes an Esri asc-file
    static void writeRasterAsASC(GeoLib::Raster const& raster, std::string const& file_name);


private:
    /// Reads the header of a Esri asc-file.
    static bool readASCHeader(std::ifstream &in, GeoLib::RasterHeader &header);

    /// Reads the header of a Surfer grd-file.
    static bool readSurferHeader(std::ifstream &in, GeoLib::RasterHeader &header,
                                 double &min, double &max);
};

/// Reads a vector of rasters given by file names. On error nothing is returned,
/// otherwise the returned vector contains pointers to the read rasters.
boost::optional<std::vector<GeoLib::Raster const*>> readRasters(
    std::vector<std::string> const& raster_paths);
} // end namespace FileIO

#endif /* ASCIIRASTERINTERFACE_H_ */
