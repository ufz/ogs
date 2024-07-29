/**
 * \file
 * \author Karsten Rink
 * \date 2014-09-10
 * \brief Implementation of the AsciiRasterInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "AsciiRasterInterface.h"

#include <fstream>
#include <tuple>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/Point.h"

namespace FileIO
{
GeoLib::Raster* AsciiRasterInterface::readRaster(std::string const& fname)
{
    std::string ext(BaseLib::getFileExtension(fname));
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);
    if (ext == ".asc")
    {
        return getRasterFromASCFile(fname);
    }
    if (ext == ".grd")
    {
        return getRasterFromSurferFile(fname);
    }
    if (ext == ".xyz")
    {
        return getRasterFromXyzFile(fname);
    }
    return nullptr;
}

/// Reads a double replacing comma by point.
static double readDoubleFromStream(std::istream& in)
{
    std::string value;
    in >> value;
    return std::strtod(BaseLib::replaceString(",", ".", value).c_str(),
                       nullptr);
}

/// Reads the header of a Esri asc-file.
/// If the return value is empty, reading was not successful.
static std::optional<GeoLib::RasterHeader> readASCHeader(std::ifstream& in)
{
    GeoLib::RasterHeader header;

    std::string tag;
    std::string value;

    in >> tag;
    if (tag == "ncols")
    {
        in >> value;
        header.n_cols = atoi(value.c_str());
    }
    else
    {
        return {};
    }

    in >> tag;
    if (tag == "nrows")
    {
        in >> value;
        header.n_rows = atoi(value.c_str());
    }
    else
    {
        return {};
    }

    header.n_depth = 1;

    in >> tag;
    if (tag == "xllcorner" || tag == "xllcenter")
    {
        header.origin[0] = readDoubleFromStream(in);
    }
    else
    {
        return {};
    }

    in >> tag;
    if (tag == "yllcorner" || tag == "yllcenter")
    {
        header.origin[1] = readDoubleFromStream(in);
    }
    else
    {
        return {};
    }
    header.origin[2] = 0;

    in >> tag;
    if (tag == "cellsize")
    {
        header.cell_size = readDoubleFromStream(in);
    }
    else
    {
        return {};
    }

    in >> tag;
    if (tag == "NODATA_value" || tag == "nodata_value")
    {
        header.no_data = readDoubleFromStream(in);
    }
    else
    {
        return {};
    }

    return header;
}

GeoLib::Raster* AsciiRasterInterface::getRasterFromASCFile(
    std::string const& fname)
{
    std::ifstream in(fname.c_str());

    if (!in.is_open())
    {
        WARN("Raster::getRasterFromASCFile(): Could not open file {:s}.",
             fname);
        return nullptr;
    }

    auto const header = readASCHeader(in);
    if (!header)
    {
        WARN(
            "Raster::getRasterFromASCFile(): Could not read header of file "
            "{:s}",
            fname);
        return nullptr;
    }

    std::vector<double> values(header->n_cols * header->n_rows);
    // read the data into the double-array
    for (std::size_t j(0); j < header->n_rows; ++j)
    {
        const std::size_t idx((header->n_rows - j - 1) * header->n_cols);
        for (std::size_t i(0); i < header->n_cols; ++i)
        {
            values[idx + i] = readDoubleFromStream(in);
        }
    }

    return new GeoLib::Raster(*header, values.begin(), values.end());
}

/// Reads the header of a Surfer grd-file with minimum and maximum values.
/// If the return value is empty, reading was not successful.
static std::optional<std::tuple<GeoLib::RasterHeader, double, double>>
readSurferHeader(std::ifstream& in)
{
    std::string tag;

    in >> tag;

    if (tag != "DSAA")
    {
        ERR("Error in readSurferHeader() - No Surfer file.");
        return {};
    }

    GeoLib::RasterHeader header;
    in >> header.n_cols >> header.n_rows;
    double min, max;
    in >> min >> max;
    header.origin[0] = min;
    header.cell_size = (max - min) / static_cast<double>(header.n_cols);

    in >> min >> max;
    header.origin[1] = min;
    header.origin[2] = 0;

    if (ceil((max - min) / static_cast<double>(header.n_rows)) ==
        ceil(header.cell_size))
    {
        header.cell_size = ceil(header.cell_size);
    }
    else
    {
        ERR("Error in readSurferHeader() - Anisotropic cellsize detected.");
        return {};
    }
    header.n_depth = 1;
    header.no_data = -9999;
    in >> min >> max;

    return {{header, min, max}};
}

GeoLib::Raster* AsciiRasterInterface::getRasterFromSurferFile(
    std::string const& fname)
{
    std::ifstream in(fname.c_str());

    if (!in.is_open())
    {
        ERR("Raster::getRasterFromSurferFile() - Could not open file {:s}",
            fname);
        return nullptr;
    }

    auto const optional_header = readSurferHeader(in);
    if (!optional_header)
    {
        ERR("Raster::getRasterFromASCFile() - could not read header of file "
            "{:s}",
            fname);
        return nullptr;
    }

    auto const [header, min, max] = *optional_header;
    std::vector<double> values(header.n_cols * header.n_rows);
    // read the data into the double-array
    for (std::size_t j(0); j < header.n_rows; ++j)
    {
        const std::size_t idx(j * header.n_cols);
        for (std::size_t i(0); i < header.n_cols; ++i)
        {
            const double val = readDoubleFromStream(in);
            values[idx + i] = (val > max || val < min) ? header.no_data : val;
        }
    }

    return new GeoLib::Raster(header, values.begin(), values.end());
}

std::optional<std::array<double, 3>> readCoordinates(std::istream& in)
{
    std::string line("");
    if (std::getline(in, line))
    {
        std::array<double, 3> coords;
        if (std::sscanf(line.c_str(), "%lf %lf %lf", &coords[0], &coords[1],
                        &coords[2]) == 3)
        {
            return coords;
        }
        else
        {
            OGS_FATAL("Could not parse coordinates in readCoordinates().");
        }
    }
    return std::nullopt;
}

GeoLib::Raster* AsciiRasterInterface::getRasterFromXyzFile(
    std::string const& fname)
{
    std::ifstream in(fname.c_str());
    if (!in.is_open())
    {
        ERR("Raster::getRasterFromXyzFile() - Could not open file {:s}", fname);
        return nullptr;
    }

    auto coords = readCoordinates(in);
    if (coords == std::nullopt)
    {
        return nullptr;
    }
    double x_min = (*coords)[0];
    double x_max = (*coords)[0];
    double y_min = (*coords)[1];
    double y_max = (*coords)[1];
    double cellsize = std::numeric_limits<double>::max();

    while ((coords = readCoordinates(in)))
    {
        double const diff = (*coords)[0] - x_min;
        if (diff > 0)
            cellsize = std::min(cellsize, diff);
        x_min = std::min((*coords)[0], x_min);
        x_max = std::max((*coords)[0], x_max);
        y_min = std::min((*coords)[1], y_min);
        y_max = std::max((*coords)[1], y_max);
    }
    in.close();

    GeoLib::RasterHeader header;
    header.cell_size = cellsize;
    header.no_data = -9999;
    header.n_cols = static_cast<std::size_t>(((x_max - x_min) / cellsize) + 1);
    header.n_rows = static_cast<std::size_t>(((y_max - y_min) / cellsize) + 1);
    header.n_depth = 1;
    header.origin[0] = x_min;
    header.origin[1] = y_min;

    std::vector<double> values(header.n_cols * header.n_rows, -9999);
    in.open(fname);
    while ((coords = readCoordinates(in)))
    {
        std::size_t idx = static_cast<std::size_t>(
            (header.n_cols * (((*coords)[1] - y_min) / cellsize)) +
            (((*coords)[0] - x_min) / cellsize));
        values[idx] = (*coords)[2];
    }
    in.close();
    return new GeoLib::Raster(header, values.begin(), values.end());
}

void AsciiRasterInterface::writeRasterAsASC(GeoLib::Raster const& raster,
                                            std::string const& file_name)
{
    GeoLib::RasterHeader header(raster.getHeader());
    MathLib::Point3d const& origin(header.origin);
    unsigned const nCols(header.n_cols);
    unsigned const nRows(header.n_rows);

    // write header
    std::ofstream out(file_name);
    out << "ncols " << nCols << "\n";
    out << "nrows " << nRows << "\n";
    auto const default_precision = out.precision();
    out.precision(std::numeric_limits<double>::digits10);
    out << "xllcorner " << origin[0] << "\n";
    out << "yllcorner " << origin[1] << "\n";
    out << "cellsize " << header.cell_size << "\n";
    out.precision(default_precision);
    out << "NODATA_value " << header.no_data << "\n";

    // write data
    for (unsigned row(0); row < nRows; ++row)
    {
        for (unsigned col(0); col < nCols - 1; ++col)
        {
            out << raster.data()[(nRows - row - 1) * nCols + col] << " ";
        }
        out << raster.data()[(nRows - row) * nCols - 1] << "\n";
    }
    out.close();
}

/// Checks if all raster files actually exist
static bool allRastersExist(std::vector<std::string> const& raster_paths)
{
    return std::all_of(raster_paths.begin(), raster_paths.end(),
                       [](std::string const& raster_path)
                       {
                           if (BaseLib::IsFileExisting(raster_path))
                           {
                               return true;
                           }
                           ERR("Opening raster file {} failed.", raster_path);
                           return false;
                       });
}

std::optional<std::vector<GeoLib::Raster const*>> readRasters(
    std::vector<std::string> const& raster_paths)
{
    if (!allRastersExist(raster_paths))
    {
        return std::nullopt;
    }

    std::vector<GeoLib::Raster const*> rasters;
    rasters.reserve(raster_paths.size());
    std::transform(raster_paths.begin(), raster_paths.end(),
                   std::back_inserter(rasters),
                   [](auto const& path)
                   { return FileIO::AsciiRasterInterface::readRaster(path); });
    return std::make_optional(rasters);
}
}  // end namespace FileIO
