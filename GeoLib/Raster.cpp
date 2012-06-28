/**
 * \file Raster.cpp
 *
 *  Created on 2011-09-07 by Thomas Fischer
 */

#include "Raster.h"

namespace GeoLib {

Raster::Raster(double cell_size, double no_data_val) :
	_cell_size(cell_size), _no_data_val(no_data_val)
{
}

void Raster::setCellSize(double cell_size)
{
	_cell_size = cell_size;
}

void Raster::setNoDataVal (double no_data_val)
{
	_no_data_val = no_data_val;
}

double* Raster::getRasterFromSurface(Surface const& sfc, size_t &n_x_pnts, size_t &n_y_pnts) const
{
	Point const& ll (sfc.getAABB().getMinPoint());
	Point const& ur (sfc.getAABB().getMaxPoint());

	n_x_pnts = static_cast<size_t>(fabs(ur[0]-ll[0]) / _cell_size)+1;
	n_y_pnts = static_cast<size_t>(fabs(ur[1]-ll[1]) / _cell_size)+1;
	const size_t n_triangles (sfc.getNTriangles());
	double *z_vals (new double[n_x_pnts*n_y_pnts]);
	if (!z_vals) {
		std::cout << "DEBUG: CreateRaster::getRaster " << n_x_pnts << " x " << n_y_pnts << " to big" << std::endl;
	}
	size_t k(0);

	for (size_t r(0); r < n_x_pnts; r++) {
		for (size_t c(0); c < n_y_pnts; c++) {
			const double test_pnt[3] = { ll[0] + r*_cell_size, ll[1] + c*_cell_size, 0};
			for (k=0; k<n_triangles; k++) {
				if (sfc[k]->containsPoint2D(test_pnt)) {
					Triangle const * const tri (sfc[k]);
					// compute coefficients c0, c1, c2 for the plane f(x,y) = c0 x + c1 y + c2
					double coeff[3] = {0.0, 0.0, 0.0};
					GeoLib::getPlaneCoefficients(*tri, coeff);
					z_vals[r*n_y_pnts+c] = coeff[0] * test_pnt[0] + coeff[1] * test_pnt[1] + coeff[2];
					break;
				}
			}
			if (k==n_triangles) {
				z_vals[r*n_y_pnts+c] = _no_data_val;
			}
		}
	}

	return z_vals;
}

Raster::~Raster()
{}

}
