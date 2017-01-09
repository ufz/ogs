/**
 * @file Raster.h
 * @author Thomas Fischer
 * @date 2011-09-07
 * @brief Definition of the GeoLib::Raster class.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Surface.h"

namespace GeoLib {

/// Contains the relevant information when handling with geoscientific raster data
struct RasterHeader
{
    std::size_t n_cols; // width
    std::size_t n_rows; // height
    MathLib::Point3d origin; // lower left corner
    double cell_size; // edge length of each pixel
    double no_data; // no data value
};

/**
 * @brief Class Raster is used for managing raster data.
 *
 * A raster consists of the meta data like number of columns and rows, the lower
 * left point, the size of a raster pixel and a value for invalid data pixels.
 * Additional the object needs the raster data itself. The raster data will be
 * copied from the constructor. The destructor will release the memory.
 */
class Raster {
public:
    typedef double const* const_iterator;
    typedef double* iterator;

    /**
     * @brief Constructor for an object of class Raster. The raster data have
     * to be handed over via input iterators. Deploying iterators has the
     * advantage that the use of the class is independent from the input
     * container.
     * @param header meta-information about the raster (height, width, etc.)
     * @param begin input iterator pointing to the first element of the data
     * @param end input iterator pointing to the last element of the data, end have to be reachable from begin
     */
    template<typename InputIterator>
    Raster(RasterHeader header, InputIterator begin, InputIterator end) :
        _header(header), _raster_data(new double[_header.n_cols*_header.n_rows])
    {
        iterator raster_it(_raster_data);
        for (InputIterator it(begin); it != end; ++it) {
            *raster_it = *it;
            raster_it++;
        }
    }

    /// Returns the complete header information
    RasterHeader const& getHeader() const { return _header; }

    /**
     * Refine the raster using scaling as a refinement parameter.
     */
    void refineRaster(std::size_t scaling);

    /**
     * Constant iterator that is pointing to the first raster pixel value.
     * @return constant iterator
     */
    const_iterator begin() const { return _raster_data; }
    /**
     * Constant iterator that is pointing to the last raster pixel value.
     * @return constant iterator
     */
    const_iterator end() const { return _raster_data + _header.n_rows*_header.n_cols; }

    /**
     * Returns the raster value at the position of the given point.
     */
    double getValueAtPoint(const MathLib::Point3d &pnt) const;

    /// Interpolates the elevation of the given point based on the 8-neighbourhood of the raster cell it is located on
    double interpolateValueAtPoint(const MathLib::Point3d &pnt) const;

    /// Checks if the given point is located within the (x,y)-extension of the raster.
    bool isPntOnRaster(MathLib::Point3d const& node) const;

    ~Raster();

    /// Creates a Raster based on a GeoLib::Surface
    static Raster* getRasterFromSurface(Surface const& sfc, double cell_size, double no_data_val = -9999);

private:
    void setCellSize(double cell_size);
    void setNoDataVal (double no_data_val);

    GeoLib::RasterHeader _header;
    double* _raster_data;
};

}
