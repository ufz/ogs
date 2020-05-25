/**
 * @file Raster.cpp
 * @author Thomas Fischer
 * @date 2011-09-07
 * @brief Implementation of the GeoLib::Raster class.
 *
 * @copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <fstream>

#include "Raster.h"

// BaseLib
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "Triangle.h"

namespace GeoLib {

void Raster::refineRaster(std::size_t scaling)
{
    auto* new_raster_data(
        new double[header_.n_rows * header_.n_cols * scaling * scaling]);

    for (std::size_t row(0); row<header_.n_rows; row++) {
        for (std::size_t col(0); col<header_.n_cols; col++) {
            const std::size_t idx(row*header_.n_cols+col);
            for (std::size_t new_row(row*scaling); new_row<(row+1)*scaling; new_row++) {
                const std::size_t idx0(new_row*header_.n_cols*scaling);
                for (std::size_t new_col(col*scaling); new_col<(col+1)*scaling; new_col++) {
                    new_raster_data[idx0+new_col] = raster_data_[idx];
                }
            }
        }
    }

    std::swap(raster_data_, new_raster_data);
    header_.cell_size /= scaling;
    header_.n_cols *= scaling;
    header_.n_rows *= scaling;

    delete [] new_raster_data;
}

Raster::~Raster()
{
    delete [] raster_data_;
}

double Raster::getValueAtPoint(const MathLib::Point3d &pnt) const
{
    if (pnt[0]>=header_.origin[0] && pnt[0]<(header_.origin[0]+(header_.cell_size*header_.n_cols)) &&
        pnt[1]>=header_.origin[1] && pnt[1]<(header_.origin[1]+(header_.cell_size*header_.n_rows)))
    {
        auto cell_x = static_cast<int>(
            std::floor((pnt[0] - header_.origin[0]) / header_.cell_size));
        auto cell_y = static_cast<int>(
            std::floor((pnt[1] - header_.origin[1]) / header_.cell_size));

        // use raster boundary values if node is outside raster due to rounding
        // errors or floating point arithmetic
        cell_x = (cell_x < 0) ? 0 : ((cell_x > static_cast<int>(header_.n_cols))
                                         ? static_cast<int>(header_.n_cols - 1)
                                         : cell_x);
        cell_y = (cell_y < 0) ? 0 : ((cell_y > static_cast<int>(header_.n_rows))
                                         ? static_cast<int>(header_.n_rows - 1)
                                         : cell_y);

        const std::size_t index = cell_y * header_.n_cols + cell_x;
        return raster_data_[index];
    }
    return header_.no_data;
}

double Raster::interpolateValueAtPoint(MathLib::Point3d const& pnt) const
{
    // position in raster
    double const xPos ((pnt[0] - header_.origin[0]) / header_.cell_size);
    double const yPos ((pnt[1] - header_.origin[1]) / header_.cell_size);
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
    std::array<double,4>  pix_val{};
    unsigned no_data_count (0);
    for (unsigned j=0; j<4; ++j)
    {
        // check if neighbour pixel is still on the raster, otherwise substitute
        // a no data value. This also allows the cast to unsigned type.
        if ((xIdx + x_nb[j]) < 0 || (yIdx + y_nb[j]) < 0 ||
            (xIdx + x_nb[j]) > (header_.n_cols - 1) ||
            (yIdx + y_nb[j]) > (header_.n_rows - 1))
        {
            pix_val[j] = header_.no_data;
        }
        else
        {
            pix_val[j] = raster_data_[static_cast<std::size_t>(yIdx + y_nb[j]) *
                                          header_.n_cols +
                                      static_cast<std::size_t>(xIdx + x_nb[j])];
        }

        // remove no data values
        if (std::fabs(pix_val[j] - header_.no_data) < std::numeric_limits<double>::epsilon())
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
            return header_.no_data;
        }

        const double norm = 1.0 / (weight[0]+weight[1]+weight[2]+weight[3]);
        std::for_each(weight.begin(), weight.end(), [&norm](double &val){val*=norm;});
    }

    // new value
    return MathLib::scalarProduct<double,4>(weight.data(), pix_val.data());
    }

bool Raster::isPntOnRaster(MathLib::Point3d const& pnt) const
{
    return !(
        (pnt[0] < header_.origin[0]) ||
        (pnt[0] > header_.origin[0] + (header_.n_cols * header_.cell_size)) ||
        (pnt[1] < header_.origin[1]) ||
        (pnt[1] > header_.origin[1] + (header_.n_rows * header_.cell_size)));
}

} // end namespace GeoLib
