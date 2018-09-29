/**
 * @file Raster.cpp
 * @author Thomas Fischer
 * @date 2011-09-07
 * @brief Implementation of the GeoLib::Raster class.
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <fstream>

#include <logog/include/logog.hpp>

#include "Raster.h"

// BaseLib
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "Triangle.h"

namespace GeoLib {

void Raster::refineRaster(std::size_t scaling)
{
    auto* new_raster_data(
        new double[_header.n_rows * _header.n_cols * scaling * scaling]);

    for (std::size_t row(0); row<_header.n_rows; row++) {
        for (std::size_t col(0); col<_header.n_cols; col++) {
            const std::size_t idx(row*_header.n_cols+col);
            for (std::size_t new_row(row*scaling); new_row<(row+1)*scaling; new_row++) {
                const std::size_t idx0(new_row*_header.n_cols*scaling);
                for (std::size_t new_col(col*scaling); new_col<(col+1)*scaling; new_col++) {
                    new_raster_data[idx0+new_col] = _raster_data[idx];
                }
            }
        }
    }

    std::swap(_raster_data, new_raster_data);
    _header.cell_size /= scaling;
    _header.n_cols *= scaling;
    _header.n_rows *= scaling;

    delete [] new_raster_data;
}

Raster::~Raster()
{
    delete [] _raster_data;
}

Raster* Raster::getRasterFromSurface(Surface const& sfc, double cell_size, double no_data_val)
{
    MathLib::Point3d const& ll(sfc.getAABB().getMinPoint());
    MathLib::Point3d const& ur(sfc.getAABB().getMaxPoint());

    const std::size_t n_cols = static_cast<std::size_t>(std::abs(ur[0]-ll[0]) / cell_size)+1;
    const std::size_t n_rows = static_cast<std::size_t>(std::abs(ur[1]-ll[1]) / cell_size)+1;
    const std::size_t n_triangles(sfc.getNumberOfTriangles());
    auto* z_vals(new double[n_cols * n_rows]);
    std::size_t k(0);

    for (std::size_t r(0); r < n_cols; r++) {
        for (std::size_t c(0); c < n_rows; c++) {
            GeoLib::Point const test_pnt = { ll[0] + r*cell_size, ll[1] + c*cell_size, 0};
            for (k=0; k<n_triangles; k++) {
                if (sfc[k]->containsPoint2D(test_pnt)) {
                    GeoLib::Triangle const * const tri (sfc[k]);
                    // compute coefficients c0, c1, c2 for the plane f(x,y) = c0 x + c1 y + c2
                    double coeff[3] = {0.0, 0.0, 0.0};
                    GeoLib::getPlaneCoefficients(*tri, coeff);
                    z_vals[r*n_rows+c] = coeff[0] * test_pnt[0] + coeff[1] * test_pnt[1] + coeff[2];
                    break;
                }
            }
            if (k==n_triangles) {
                z_vals[r*n_rows+c] = no_data_val;
            }
        }
    }

    RasterHeader header = { std::size_t(n_cols),  std::size_t(n_rows), 1,
        MathLib::Point3d(ll), cell_size, static_cast<double>(-9999) };
    return new Raster(header, z_vals, z_vals+n_cols*n_rows);
}

double Raster::getValueAtPoint(const MathLib::Point3d &pnt) const
{
    if (pnt[0]>=_header.origin[0] && pnt[0]<(_header.origin[0]+(_header.cell_size*_header.n_cols)) &&
        pnt[1]>=_header.origin[1] && pnt[1]<(_header.origin[1]+(_header.cell_size*_header.n_rows)))
    {
        auto cell_x = static_cast<int>(
            std::floor((pnt[0] - _header.origin[0]) / _header.cell_size));
        auto cell_y = static_cast<int>(
            std::floor((pnt[1] - _header.origin[1]) / _header.cell_size));

        // use raster boundary values if node is outside raster due to rounding
        // errors or floating point arithmetic
        cell_x = (cell_x < 0) ? 0 : ((cell_x > static_cast<int>(_header.n_cols))
                                         ? static_cast<int>(_header.n_cols - 1)
                                         : cell_x);
        cell_y = (cell_y < 0) ? 0 : ((cell_y > static_cast<int>(_header.n_rows))
                                         ? static_cast<int>(_header.n_rows - 1)
                                         : cell_y);

        const std::size_t index = cell_y * _header.n_cols + cell_x;
        return _raster_data[index];
    }
    return _header.no_data;
}

double Raster::interpolateValueAtPoint(MathLib::Point3d const& pnt) const
{
    // position in raster
    double const xPos ((pnt[0] - _header.origin[0]) / _header.cell_size);
    double const yPos ((pnt[1] - _header.origin[1]) / _header.cell_size);
    // raster cell index
    double const xIdx (std::floor(xPos));    //carry out computions in double
    double const yIdx (std::floor(yPos));    //  so not to over- or underflow.

    // weights for bilinear interpolation
    double const xShift = std::fabs((xPos - xIdx) - 0.5);
    double const yShift = std::fabs((yPos - yIdx) - 0.5);
    std::array<double,4> weight = {{ (1-xShift)*(1-yShift), xShift*(1-yShift), xShift*yShift, (1-xShift)*yShift }};

    // neighbors to include in interpolation
    int const xShiftIdx = (xPos - xIdx >= 0.5) ? 1 : -1;
    int const yShiftIdx = (yPos - yIdx >= 0.5) ? 1 : -1;
    std::array<int,4> const x_nb = {{ 0, xShiftIdx, xShiftIdx, 0 }};
    std::array<int,4> const y_nb = {{ 0, 0, yShiftIdx, yShiftIdx }};

    // get pixel values
    std::array<double,4>  pix_val;
    unsigned no_data_count (0);
    for (unsigned j=0; j<4; ++j)
    {
        // check if neighbour pixel is still on the raster, otherwise substitute
        // a no data value. This also allows the cast to unsigned type.
        if ( (xIdx + x_nb[j]) < 0 ||
             (yIdx + y_nb[j]) < 0 ||
             (xIdx + x_nb[j]) > (_header.n_cols-1) ||
             (yIdx + y_nb[j]) > (_header.n_rows-1) )
            pix_val[j] = _header.no_data;
        else
            pix_val[j] = _raster_data[
                static_cast<std::size_t>(yIdx + y_nb[j]) * _header.n_cols +
                static_cast<std::size_t>(xIdx + x_nb[j])];

        // remove no data values
        if (std::fabs(pix_val[j] - _header.no_data) < std::numeric_limits<double>::epsilon())
        {
            weight[j] = 0;
            no_data_count++;
        }
    }

    // adjust weights if necessary
    if (no_data_count > 0)
    {
        if (no_data_count == 4) // if there is absolutely no data just use the default value
            return _header.no_data;

        const double norm = 1.0 / (weight[0]+weight[1]+weight[2]+weight[3]);
        std::for_each(weight.begin(), weight.end(), [&norm](double &val){val*=norm;});
    }

    // new value
    return MathLib::scalarProduct<double,4>(weight.data(), pix_val.data());
    }

bool Raster::isPntOnRaster(MathLib::Point3d const& pnt) const
{
    return !(
        (pnt[0] < _header.origin[0]) ||
        (pnt[0] > _header.origin[0] + (_header.n_cols * _header.cell_size)) ||
        (pnt[1] < _header.origin[1]) ||
        (pnt[1] > _header.origin[1] + (_header.n_rows * _header.cell_size)));
}

} // end namespace GeoLib
