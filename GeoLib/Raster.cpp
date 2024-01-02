/**
 * \file
 * \author Thomas Fischer
 * \date 2011-09-07
 * \brief Implementation of the GeoLib::Raster class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Raster.h"

#include <fstream>

// BaseLib
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"
#include "Triangle.h"

namespace GeoLib
{
void Raster::refineRaster(std::size_t scaling)
{
    std::vector<double> new_raster_data(_header.n_rows * _header.n_cols *
                                        scaling * scaling);

    for (std::size_t row(0); row < _header.n_rows; row++)
    {
        for (std::size_t col(0); col < _header.n_cols; col++)
        {
            const std::size_t idx(row * _header.n_cols + col);
            for (std::size_t new_row(row * scaling);
                 new_row < (row + 1) * scaling;
                 new_row++)
            {
                const std::size_t idx0(new_row * _header.n_cols * scaling);
                for (std::size_t new_col(col * scaling);
                     new_col < (col + 1) * scaling;
                     new_col++)
                {
                    new_raster_data[idx0 + new_col] = _raster_data[idx];
                }
            }
        }
    }

    std::swap(_raster_data, new_raster_data);
    _header.cell_size /= scaling;
    _header.n_cols *= scaling;
    _header.n_rows *= scaling;
}

double Raster::getValueAtPoint(const MathLib::Point3d& pnt) const
{
    if (pnt[0] >= _header.origin[0] &&
        pnt[0] < (_header.origin[0] + (_header.cell_size * _header.n_cols)) &&
        pnt[1] >= _header.origin[1] &&
        pnt[1] < (_header.origin[1] + (_header.cell_size * _header.n_rows)))
    {
        auto cell_x = static_cast<int>(
            std::floor((pnt[0] - _header.origin[0]) / _header.cell_size));
        auto cell_y = static_cast<int>(
            std::floor((pnt[1] - _header.origin[1]) / _header.cell_size));

        // use raster boundary values if node is outside raster due to rounding
        // errors or floating point arithmetic
        cell_x = (cell_x < 0) ? 0
                              : ((cell_x > static_cast<int>(_header.n_cols))
                                     ? static_cast<int>(_header.n_cols - 1)
                                     : cell_x);
        cell_y = (cell_y < 0) ? 0
                              : ((cell_y > static_cast<int>(_header.n_rows))
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
    double const xPos((pnt[0] - _header.origin[0]) / _header.cell_size);
    double const yPos((pnt[1] - _header.origin[1]) / _header.cell_size);
    // raster cell index
    double const xIdx(std::floor(xPos));  // carry out computions in double
    double const yIdx(std::floor(yPos));  //  so not to over- or underflow.

    // weights for bilinear interpolation
    double const xShift = std::abs((xPos - xIdx) - 0.5);
    double const yShift = std::abs((yPos - yIdx) - 0.5);
    Eigen::Vector4d weight = {(1 - xShift) * (1 - yShift),
                              xShift * (1 - yShift), xShift * yShift,
                              (1 - xShift) * yShift};

    // neighbors to include in interpolation
    int const xShiftIdx = (xPos - xIdx >= 0.5) ? 1 : -1;
    int const yShiftIdx = (yPos - yIdx >= 0.5) ? 1 : -1;
    std::array<int, 4> const x_nb = {{0, xShiftIdx, xShiftIdx, 0}};
    std::array<int, 4> const y_nb = {{0, 0, yShiftIdx, yShiftIdx}};

    // get pixel values
    Eigen::Vector4d pix_val{};
    unsigned no_data_count(0);
    for (unsigned j = 0; j < 4; ++j)
    {
        // check if neighbour pixel is still on the raster, otherwise substitute
        // a no data value. This also allows the cast to unsigned type.
        if ((xIdx + x_nb[j]) < 0 || (yIdx + y_nb[j]) < 0 ||
            (xIdx + x_nb[j]) > (_header.n_cols - 1) ||
            (yIdx + y_nb[j]) > (_header.n_rows - 1))
        {
            pix_val[j] = _header.no_data;
        }
        else
        {
            pix_val[j] = _raster_data[static_cast<std::size_t>(yIdx + y_nb[j]) *
                                          _header.n_cols +
                                      static_cast<std::size_t>(xIdx + x_nb[j])];
        }

        // remove no data values
        if (std::abs(pix_val[j] - _header.no_data) <
            std::numeric_limits<double>::epsilon())
        {
            weight[j] = 0;
            no_data_count++;
        }
    }

    // adjust weights if necessary
    if (no_data_count > 0)
    {
        if (no_data_count == 4)
        {  // if there is absolutely no data just use the default value
            return _header.no_data;
        }

        weight /= weight.sum();
    }

    // new value
    return weight.dot(pix_val);
}

bool Raster::isPntOnRaster(MathLib::Point3d const& pnt) const
{
    return !(
        (pnt[0] < _header.origin[0]) ||
        (pnt[0] > _header.origin[0] + (_header.n_cols * _header.cell_size)) ||
        (pnt[1] < _header.origin[1]) ||
        (pnt[1] > _header.origin[1] + (_header.n_rows * _header.cell_size)));
}

}  // end namespace GeoLib
