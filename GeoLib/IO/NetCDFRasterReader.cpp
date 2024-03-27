/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GeoLib/IO/NetCDFRasterReader.h"

#include <filesystem>
#include <numeric>
#ifdef OGS_USE_NETCDF
#include <netcdf>
#endif

#include "BaseLib/ConfigTree.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/AABB.h"
#include "GeoLib/IO/AsciiRasterInterface.h"

namespace
{
#ifdef OGS_USE_NETCDF

GeoLib::RasterHeader readHeaderFromNetCDF(
    std::multimap<std::string, netCDF::NcVar> const& variables)
{
    GeoLib::RasterHeader header;

    // remark: netcdf raster starts in the left upper corner
    double cell_size_y = 0;
    // first read:
    // - origin, cell_size from GeoTransform attribute
    // - header.no_data from _FillValue attribute
    for (auto [variable_name, variable] : variables)
    {
        if (variable.isNull())
        {
            OGS_FATAL("Variable '{}' not found in file.", variable_name);
        }
        auto const& attributes = variable.getAtts();
        for (auto [name, attribute] : attributes)
        {
            if (name == "GeoTransform")
            {
                std::vector<char> attribute_values;
                attribute_values.resize(attribute.getAttLength());
                attribute.getValues(attribute_values.data());
                std::string const str_attribute_values(attribute_values.begin(),
                                                       attribute_values.end());
                auto const strings = BaseLib::splitString(str_attribute_values);
                header.origin[0] = std::stod(strings[0]);  // origin x
                header.cell_size = std::stod(strings[1]);  // cell size x
                header.origin[1] = std::stod(strings[3]);  // origin y
                cell_size_y = std::stod(strings[5]);       // cell size y
                DBUG(
                    "GeoTransform attribute values: origin_x: {}, origin_y: "
                    "{}, cell_size_x: {}, cell_size_y: {}",
                    header.origin[0], header.origin[1], header.cell_size,
                    cell_size_y);
            }
            if (name == "_FillValue")
            {
                attribute.getValues(&header.no_data);
            }
        }
    }

    return header;
}

std::vector<double> readDataFromNetCDF(
    std::multimap<std::string, netCDF::NcVar> const& variables,
    std::string_view const var_name, std::size_t const dimension_number,
    GeoLib::MinMaxPoints const& min_max_points, GeoLib::RasterHeader& header)
{
    // second read raster header values
    std::vector<double> raster_data;
    for (auto [variable_name, variable] : variables)
    {
        // uncomment for debugging
        // DBUG("checking variable '{}'", variable_name);
        if (variable_name == var_name)
        {
            auto const n_dims = variable.getDimCount();
            if (n_dims == 0)
            {
                continue;
            }
            // DBUG("variable '{}' has {} dimensions", variable_name, n_dims);

            std::vector<std::size_t> counts(n_dims, 1);
            for (int i = 0; i < n_dims; ++i)
            {
                counts[i] = variable.getDim(i).getSize();
            }
            // DBUG("counts {}", fmt::join(counts, ", "));
            std::vector<std::size_t> start(n_dims, 0);
            if (dimension_number > counts[0])
            {
                OGS_FATAL(
                    "variable '{}' has {} dimensions, requested dimension "
                    "number is {}",
                    variable_name, counts[0], dimension_number);
            }
            // time dimension
            start[0] = dimension_number;  // number of time step
            counts[0] = 1;                // read one time step

            // y dimension - netcdf dimension number 1
            // With the counts[1] and cell_size the 'usual' lower y-coordinate
            // of the origin is computed.
            // DBUG("reset header.origin[1]: original y0: {}, new y0: {} ",
            //     header.origin[1],
            //     header.origin[1] - counts[1] * header.cell_size);
            header.origin[1] -= counts[1] * header.cell_size;

            counts[2] = static_cast<int>(
                std::ceil((min_max_points.max[0] - min_max_points.min[0]) /
                          header.cell_size) +
                1);
            counts[1] = static_cast<int>(
                std::ceil((min_max_points.max[1] - min_max_points.min[1]) /
                          header.cell_size) +
                1);
            // x dimension - netcdf dimension number 2
            start[2] = static_cast<int>(std::floor(
                (min_max_points.min[0] - header.origin[0]) / header.cell_size));
            // y dimension - netcdf dimension number 1
            start[1] = static_cast<int>(std::floor(
                (min_max_points.min[1] - header.origin[1]) / header.cell_size));

            // DBUG("cut-out raster: pixel in x dimension: {} ", counts[2]);
            // DBUG("cut-out raster: pixel in y dimension: {} ", counts[1]);
            // DBUG("cut-out raster: start index x dimension: {}", start[2]);
            // DBUG("cut-out raster: start index y dimension: {}", start[1]);

            header.n_cols = counts[2];
            header.n_rows = counts[1];
            header.origin[0] += start[2] * header.cell_size;
            // reset header y origin for cut out raster
            header.origin[1] += start[1] * header.cell_size;
            std::size_t const total_length =
                std::accumulate(counts.cbegin(), counts.cend(), 1,
                                std::multiplies<std::size_t>());
            raster_data.resize(total_length);
            variable.getVar(start, counts, raster_data.data());

            std::replace(raster_data.begin(), raster_data.end(), header.no_data,
                         0.0);
        }
    }
    return raster_data;
}

std::unique_ptr<GeoLib::Raster> readNetCDF(
    std::filesystem::path const& filepath,
    std::string_view const var_name,
    std::size_t const dimension_number,
    GeoLib::MinMaxPoints const& min_max_points)
{
    netCDF::NcFile dataset(filepath.string(), netCDF::NcFile::read);
    if (dataset.isNull())
    {
        OGS_FATAL("Error opening file '{}'.", filepath.string());
    }

    auto const& variables = dataset.getVars();
    GeoLib::RasterHeader header = readHeaderFromNetCDF(variables);
    std::vector<double> raster_data = readDataFromNetCDF(
        variables, var_name, dimension_number, min_max_points, header);

    return std::make_unique<GeoLib::Raster>(header, raster_data.begin(),
                                            raster_data.end());
}
#endif

GeoLib::NamedRaster readRasterFromFile(
    std::filesystem::path const& path,
    std::filesystem::path filename,
    std::string const& var_name,
    std::size_t const dimension_number,
    [[maybe_unused]] GeoLib::MinMaxPoints const& min_max_points)
{
    INFO("readRasterFromFile: '{}/{}'", path.string(), filename.string());

    if (filename.extension() == ".nc")
    {
#ifdef OGS_USE_NETCDF
        auto raster = readNetCDF(path / filename, var_name, dimension_number,
                                 min_max_points);

        return GeoLib::NamedRaster{filename.replace_extension().string() + "_" +
                                       var_name + "_" +
                                       std::to_string(dimension_number),
                                   std::move(raster)};
#else
        OGS_FATAL("OGS was not build with NetCDF support. Can not read {}",
                  (path / filename).string());
#endif
    }
    auto raster = std::unique_ptr<GeoLib::Raster>(
        FileIO::AsciiRasterInterface::readRaster((path / filename).string()));
    if (raster == nullptr)
    {
        OGS_FATAL("Could not read raster from file '{}'.",
                  (path / filename).string());
    }
    return GeoLib::NamedRaster{filename.replace_extension().string() + "_" +
                                   var_name + "_" +
                                   std::to_string(dimension_number),
                               std::move(raster)};
}
}  // anonymous namespace

namespace GeoLib::IO
{
GeoLib::NamedRaster readRaster(BaseLib::ConfigTree const& raster_config,
                               std::string const& raster_directory,
                               GeoLib::MinMaxPoints const& min_max_points)
{
    auto const file_name =
        //! \ogs_file_param{prj__rasters__raster__file}
        raster_config.getConfigParameter<std::string>("file");
    auto const variable_name =
        //! \ogs_file_param{prj__rasters__raster__variable}
        raster_config.getConfigParameter<std::string>("variable");
    auto const dimension =
        //! \ogs_file_param{prj__rasters__raster__dimension}
        raster_config.getConfigParameter<std::size_t>("dimension", 1);
    return readRasterFromFile(raster_directory, file_name, variable_name,
                              dimension, min_max_points);
}
}  // namespace GeoLib::IO
