/**
 * \file
 * \author Thomas Fischer
 * \date 2011-09-07
 * \brief Definition of the GeoLib::Raster class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
#include <stdexcept>
#include <utility>

#include "BaseLib/Error.h"
#include "MathLib/Point3d.h"

namespace GeoLib
{

/// Contains the relevant information when storing a geoscientific raster data
struct RasterHeader final
{
    std::size_t n_cols;       // width
    std::size_t n_rows;       // height
    std::size_t n_depth;      // depth (for 3d image)
    MathLib::Point3d origin;  // lower left corner
    double cell_size;         // edge length of each pixel
    double no_data;           // no data value

    bool operator==(RasterHeader const&) const = default;
};

/**
 * \brief Class Raster is used for managing raster data.
 *
 * A raster consists of the meta data like number of columns and rows, the
 * lower left point, the size of a raster pixel and a value for invalid data
 * pixels. Additional the object needs the raster data itself. The raster
 * data will be copied from the constructor. The destructor will release the
 * memory.
 */
class Raster final
{
public:
    /**
     * \brief Constructor for an object of class Raster. The raster data
     * have to be handed over via input iterators. Deploying iterators has
     * the advantage that the use of the class is independent from the input
     * container.
     * \param header meta-information about the raster (height, width, etc.)
     * \param begin input iterator pointing to the first element of the data
     * \param end input iterator pointing to the last element of the data,
     * end have to be reachable from begin
     */
    template <typename InputIterator>
    Raster(RasterHeader header, InputIterator begin, InputIterator end)
        : _header(std::move(header)),
          _raster_data(_header.n_cols * _header.n_rows)
    {
        unsigned long const number_of_input_values =
            static_cast<unsigned long>(std::distance(begin, end));
        if (number_of_input_values != _header.n_cols * _header.n_rows)
        {
            throw std::out_of_range(
                "Number of raster data mismatch, need " +
                std::to_string(_header.n_cols * _header.n_rows) +
                " values, but got " + std::to_string(number_of_input_values));
        }
        std::copy(begin, end, _raster_data.begin());
    }

    Raster(Raster const&) = default;
    Raster(Raster&&) = default;
    Raster& operator=(Raster const&) = default;
    Raster& operator=(Raster&&) = default;
    ~Raster() = default;

    /// Returns the complete header information
    RasterHeader const& getHeader() const { return _header; }

    /**
     * Refine the raster using scaling as a refinement parameter.
     */
    void refineRaster(std::size_t scaling);

    /**
     * Constant iterator that is pointing to the first raster pixel value.
     * \return constant iterator
     */
    std::vector<double>::const_iterator begin() const
    {
        return _raster_data.begin();
    }

    /**
     * Iterator to the element following the last pixel element of the
     * raster
     * \return constant iterator
     */
    std::vector<double>::const_iterator end() const
    {
        return _raster_data.end();
    }

    double const* data() const { return _raster_data.data(); }

    /**
     * Access the pixel specified by row, col.
     */
    double const& operator()(std::size_t const row, std::size_t const col) const
    {
        if (row >= _header.n_rows || col >= _header.n_cols)
        {
            OGS_FATAL(
                "Raster pixel ({}, {}) doesn't exist. Raster size is {} x {}.",
                row, col, _header.n_rows, _header.n_cols);
        }
        return _raster_data[(_header.n_rows - 1 - row) * _header.n_cols + col];
    }
    double& operator()(std::size_t const row, std::size_t const col)
    {
        return const_cast<double&>(std::as_const(*this)(row, col));
    }

    /**
     * Returns the raster value at the position of the given point.
     */
    double getValueAtPoint(const MathLib::Point3d& pnt) const;

    /// Interpolates the elevation of the given point based on the
    /// 8-neighbourhood of the raster cell it is located on
    double interpolateValueAtPoint(const MathLib::Point3d& pnt) const;

    /// Checks if the given point is located within the (x,y)-extension of
    /// the raster.
    bool isPntOnRaster(MathLib::Point3d const& pnt) const;

    bool operator==(Raster const&) const = default;

private:
    void setCellSize(double cell_size);
    void setNoDataVal(double no_data_val);

    GeoLib::RasterHeader _header;
    std::vector<double> _raster_data;
};

struct NamedRaster
{
    std::string raster_name;
    std::unique_ptr<GeoLib::Raster> raster;
};

}  // namespace GeoLib
