/**
 * @file AsciiRasterInterface.h
 * @author Karsten Rink
 * @date 2014-09-10
 * @brief Definition of the AsciiRasterInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef ASCIIRASTERINTERFACE_H_
#define ASCIIRASTERINTERFACE_H_

#include <fstream>

namespace GeoLib {
    class Raster;
}

namespace FileIO {

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
    static bool readASCHeader(std::ifstream &in, std::size_t &n_cols, std::size_t &n_rows,
                              double &xllcorner, double &yllcorner, double &cell_size, double &no_data_val);

    /// Reads the header of a Surfer grd-file.
    static bool readSurferHeader(std::ifstream &in, size_t &n_cols, std::size_t &n_rows,
                                 double &xllcorner, double &yllcorner, double &cell_size, double &min, double &max);
};

}

#endif /* ASCIIRASTERINTERFACE_H_ */
