/**
 * @file AsciiRasterInterface.cpp
 * @author Karsten Rink
 * @date 2014-09-10
 * @brief Implementation of the AsciiRasterInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "AsciiRasterInterface.h"

#include <logog/include/logog.hpp>
#include <boost/optional.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/Raster.h"

namespace FileIO
{

GeoLib::Raster* AsciiRasterInterface::readRaster(std::string const& fname)
{
    std::string ext (BaseLib::getFileExtension(fname));
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);
    if (ext == "asc")
        return getRasterFromASCFile(fname);
    if (ext == "grd")
        return getRasterFromSurferFile(fname);
    return nullptr;
}

GeoLib::Raster* AsciiRasterInterface::getRasterFromASCFile(std::string const& fname)
{
    std::ifstream in(fname.c_str());

    if (!in.is_open()) {
        WARN("Raster::getRasterFromASCFile(): Could not open file %s.", fname.c_str());
        return nullptr;
    }

    // header information
    GeoLib::RasterHeader header;
    if (readASCHeader(in, header)) {
        auto* values = new double[header.n_cols * header.n_rows];
        std::string s;
        // read the data into the double-array
        for (std::size_t j(0); j < header.n_rows; ++j) {
            const std::size_t idx ((header.n_rows - j - 1) * header.n_cols);
            for (std::size_t i(0); i < header.n_cols; ++i) {
                in >> s;
                values[idx + i] = strtod(
                    BaseLib::replaceString(",", ".", s).c_str(), nullptr);
            }
        }
        in.close();
        GeoLib::Raster *raster(new GeoLib::Raster(header, values, values+header.n_cols*header.n_rows));
        delete [] values;
        return raster;
    }
    WARN("Raster::getRasterFromASCFile(): Could not read header of file %s",
         fname.c_str());
    return nullptr;
}

bool AsciiRasterInterface::readASCHeader(std::ifstream &in, GeoLib::RasterHeader &header)
{
    std::string tag, value;

    in >> tag;
    if (tag == "ncols")
    {
        in >> value;
        header.n_cols = atoi(value.c_str());
    } else return false;

    in >> tag;
    if (tag == "nrows")
    {
        in >> value;
        header.n_rows = atoi(value.c_str());
    } else return false;

    header.n_depth = 1;

    in >> tag;
    if (tag == "xllcorner" || tag == "xllcenter")
    {
        in >> value;
        header.origin[0] =
            strtod(BaseLib::replaceString(",", ".", value).c_str(), nullptr);
    } else return false;

    in >> tag;
    if (tag == "yllcorner" || tag == "yllcenter")
    {
        in >> value;
        header.origin[1] =
            strtod(BaseLib::replaceString(",", ".", value).c_str(), nullptr);
    } else return false;
    header.origin[2] = 0;

    in >> tag;
    if (tag == "cellsize")
    {
        in >> value;
        header.cell_size =
            strtod(BaseLib::replaceString(",", ".", value).c_str(), nullptr);
    } else return false;

    in >> tag;
    if (tag == "NODATA_value" || tag == "nodata_value")
    {
        in >> value;
        header.no_data =
            strtod(BaseLib::replaceString(",", ".", value).c_str(), nullptr);
    } else return false;

    return true;
}

GeoLib::Raster* AsciiRasterInterface::getRasterFromSurferFile(std::string const& fname)
{
    std::ifstream in(fname.c_str());

    if (!in.is_open()) {
        ERR("Raster::getRasterFromSurferFile() - Could not open file %s", fname.c_str());
        return nullptr;
    }

    // header information
    GeoLib::RasterHeader header;
    double min(0.0), max(0.0);

    if (readSurferHeader(in, header, min, max))
    {
        const double no_data_val (min-1);
        auto* values = new double[header.n_cols * header.n_rows];
        std::string s;
        // read the data into the double-array
        for (std::size_t j(0); j < header.n_rows; ++j)
        {
            const std::size_t idx (j * header.n_cols);
            for (std::size_t i(0); i < header.n_cols; ++i)
            {
                in >> s;
                const double val(strtod(
                    BaseLib::replaceString(",", ".", s).c_str(), nullptr));
                values[idx+i] = (val > max || val < min) ? no_data_val : val;
            }
        }
        in.close();
        GeoLib::Raster *raster(new GeoLib::Raster(header, values, values+header.n_cols*header.n_rows));
        delete [] values;
        return raster;
    }
    ERR("Raster::getRasterFromASCFile() - could not read header of file %s",
        fname.c_str());
    return nullptr;
}

bool AsciiRasterInterface::readSurferHeader(
    std::ifstream &in, GeoLib::RasterHeader &header, double &min, double &max)
{
    std::string tag;

    in >> tag;

    if (tag != "DSAA")
    {
        ERR("Error in readSurferHeader() - No Surfer file.");
        return false;
    }

    in >> header.n_cols >> header.n_rows;
    in >> min >> max;
    header.origin[0] = min;
    header.cell_size = (max - min) / static_cast<double>(header.n_cols);

    in >> min >> max;
    header.origin[1] = min;
    header.origin[2] = 0;

    if (ceil((max - min) / static_cast<double>(header.n_rows)) ==
        ceil(header.cell_size))
        header.cell_size = ceil(header.cell_size);
    else
    {
        ERR("Error in readSurferHeader() - Anisotropic cellsize detected.");
        return false;
    }
    header.n_depth = 1;
    header.no_data = min - 1;
    in >> min >> max;

    return true;
}

void AsciiRasterInterface::writeRasterAsASC(GeoLib::Raster const& raster, std::string const& file_name)
{
    GeoLib::RasterHeader header (raster.getHeader());
    MathLib::Point3d const& origin (header.origin);
    unsigned const nCols (header.n_cols);
    unsigned const nRows (header.n_rows);

    // write header
    std::ofstream out(file_name);
    out << "ncols " << nCols << "\n";
    out << "nrows " << nRows << "\n";
    out << "xllcorner " << origin[0] << "\n";
    out << "yllcorner " << origin[1] << "\n";
    out << "cellsize " <<  header.cell_size << "\n";
    out << "NODATA_value " << header.no_data << "\n";

    // write data
    double const*const elevation(raster.begin());
    for (unsigned row(0); row < nRows; ++row)
    {
        for (unsigned col(0); col < nCols; ++col)
        {
            out << elevation[(nRows-row-1) * nCols + col] << " ";
        }
        out << "\n";
    }
    out.close();
}


/// Checks if all raster files actually exist
static bool allRastersExist(std::vector<std::string> const& raster_paths)
{
    for (const auto& raster_path : raster_paths)
    {
        std::ifstream file_stream(raster_path, std::ifstream::in);
        if (!file_stream.good()) return false;
        file_stream.close();
    }
    return true;
}

boost::optional<std::vector<GeoLib::Raster const*>> readRasters(
    std::vector<std::string> const& raster_paths)
{
    if (!allRastersExist(raster_paths)) return boost::none;

    std::vector<GeoLib::Raster const*> rasters;
    rasters.reserve(raster_paths.size());
    for (auto const& path : raster_paths)
        rasters.push_back(FileIO::AsciiRasterInterface::readRaster(path));
    return boost::make_optional(rasters);
}
} // end namespace FileIO
