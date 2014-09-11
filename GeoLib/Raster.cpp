/**
 * @file Raster.cpp
 * @author Thomas Fischer
 * @date 2011-09-07
 * @brief Implementation of the GeoLib::Raster class.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <fstream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "Raster.h"

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

namespace GeoLib {

void Raster::refineRaster(std::size_t scaling)
{
	double *new_raster_data(new double[_n_rows*_n_cols*scaling*scaling]);

	for (std::size_t row(0); row<_n_rows; row++) {
		for (std::size_t col(0); col<_n_cols; col++) {
			const size_t idx(row*_n_cols+col);
			for (std::size_t new_row(row*scaling); new_row<(row+1)*scaling; new_row++) {
				const size_t idx0(new_row*_n_cols*scaling);
				for (std::size_t new_col(col*scaling); new_col<(col+1)*scaling; new_col++) {
					new_raster_data[idx0+new_col] = _raster_data[idx];
				}
			}
		}
	}

	std::swap(_raster_data, new_raster_data);
	_cell_size /= scaling;
	_n_cols *= scaling;
	_n_rows *= scaling;

	delete [] new_raster_data;
}

Raster::~Raster()
{
	delete [] _raster_data;
}

void Raster::setCellSize(double cell_size)
{
	_cell_size = cell_size;
}

void Raster::setNoDataVal (double no_data_val)
{
	_no_data_val = no_data_val;
}

GeoLib::Point const& Raster::getOrigin() const
{
	return _ll_pnt;
}

Raster* Raster::getRasterFromSurface(Surface const& sfc, double cell_size, double no_data_val)
{
	Point const& ll (sfc.getAABB().getMinPoint());
	Point const& ur (sfc.getAABB().getMaxPoint());

	const std::size_t n_cols = static_cast<size_t>(fabs(ur[0]-ll[0]) / cell_size)+1;
	const std::size_t n_rows = static_cast<size_t>(fabs(ur[1]-ll[1]) / cell_size)+1;
	const size_t n_triangles (sfc.getNTriangles());
	double *z_vals (new double[n_cols*n_rows]);
	size_t k(0);

	for (size_t r(0); r < n_cols; r++) {
		for (size_t c(0); c < n_rows; c++) {
			const double test_pnt[3] = { ll[0] + r*cell_size, ll[1] + c*cell_size, 0};
			for (k=0; k<n_triangles; k++) {
				if (sfc[k]->containsPoint2D(test_pnt)) {
					Triangle const * const tri (sfc[k]);
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

	return new Raster(n_cols, n_rows, ll[0], ll[1], cell_size, z_vals, z_vals+n_cols*n_rows ,-9999);
}

double Raster::getValueAtPoint(const GeoLib::Point &pnt) const
{
	if (pnt[0]>=_ll_pnt[0] && pnt[0]<(_ll_pnt[0]+(_cell_size*_n_cols)) && 
		pnt[1]>=_ll_pnt[1] && pnt[1]<(_ll_pnt[1]+(_cell_size*_n_rows)))
	{
		int cell_x = static_cast<int>(floor((pnt[0] - _ll_pnt[0])/_cell_size));
		int cell_y = static_cast<int>(floor((pnt[1] - _ll_pnt[1])/_cell_size));

		// use raster boundary values if node is outside raster due to rounding errors or floating point arithmetic
		cell_x = (cell_x < 0) ?  0 : ((cell_x > static_cast<int>(_n_cols)) ? (_n_cols-1) : cell_x);
		cell_y = (cell_y < 0) ?  0 : ((cell_y > static_cast<int>(_n_rows)) ? (_n_rows-1) : cell_y);

		const std::size_t index = cell_y*_n_cols+cell_x;
		return _raster_data[index];
	}
	return _no_data_val;
}

double Raster::interpolateValueAtPoint(GeoLib::Point const& pnt) const
{
    // position in raster
    double const xPos ((pnt[0] - _ll_pnt[0]) / _cell_size);
    double const yPos ((pnt[1] - _ll_pnt[1]) / _cell_size);
    // raster cell index
    std::size_t const xIdx (static_cast<size_t>(floor(xPos)));
    std::size_t const yIdx (static_cast<size_t>(floor(yPos)));

    // weights for bilinear interpolation
    double const half_delta = 0.5*_cell_size;
    double const xShift = fabs(xPos-(xIdx+half_delta)) / _cell_size;
    double const yShift = fabs(yPos-(yIdx+half_delta)) / _cell_size;
    std::array<double,4> weight = {{ (1-xShift)*(1-xShift), xShift*(1-yShift), xShift*yShift, (1-xShift)*yShift }};

    // neightbors to include in interpolation
    int const xShiftIdx = (xPos-xIdx-half_delta>=0) ? 1 : -1;
    int const yShiftIdx = (yPos-yIdx-half_delta>=0) ? 1 : -1;
    std::array<int,4> const x_nb = {{ 0, xShiftIdx, xShiftIdx, 0 }};
    std::array<int,4> const y_nb = {{ 0, 0, yShiftIdx, yShiftIdx }};

    // get pixel values
    std::array<double,4>  pix_val;
    unsigned no_data_count (0);
    for (unsigned j=0; j<4; ++j)
    {
        pix_val[j] = _raster_data[(yIdx + y_nb[j])*_n_cols + (xIdx + x_nb[j])];
        if (fabs(pix_val[j] - _no_data_val) < std::numeric_limits<double>::epsilon())
        {
            weight[j] = 0;
            no_data_count++;
        }
    }

    // adjust weights if necessary
    if (no_data_count > 0)
    {
        if (no_data_count == 4) // if there is absolutely no data just use the default value
            return _no_data_val;

        const double norm = (double)(4)/(4-no_data_count);
        std::for_each(weight.begin(), weight.end(), [&norm](double &val){val*=norm;});
    }

    // new value
    return MathLib::scalarProduct<double,4>(weight.data(), pix_val.data());
}

bool Raster::isPntOnRaster(GeoLib::Point const& pnt) const
{
    if ((pnt[0]<_ll_pnt[0]) || (pnt[0]>_ll_pnt[0]+(_n_cols*_cell_size)) || 
        (pnt[1]<_ll_pnt[1]) || (pnt[1]>_ll_pnt[1]+(_n_rows*_cell_size)))
		return false;
    return true;
}

} // end namespace GeoLib
