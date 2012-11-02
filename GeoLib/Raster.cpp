/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Raster.cpp
 *
 * Created on 2011-09-07 by Thomas Fischer
 */

#include <fstream>

#include "Raster.h"

// BaseLib
#include "StringTools.h"

namespace GeoLib {

Raster::Raster(std::size_t n_cols, std::size_t n_rows, double xllcorner, double yllcorner,
				double cell_size, double no_data_val, double* raster_data) :
	_n_cols(n_cols), _n_rows(n_rows), _ll_pnt(xllcorner, yllcorner, 0.0),
	_cell_size(cell_size), _no_data_val(no_data_val), _raster_data(raster_data)
{}

void Raster::refineRaster(std::size_t n_cols, std::size_t n_rows)
{
	if (n_rows <= _n_rows || n_cols <= _n_cols) return;

	std::size_t row_blk_size(n_rows / _n_rows);
	std::size_t col_blk_size(n_cols / _n_cols);
	double *new_raster_data(new double[n_rows*n_cols]);

	for (std::size_t row(0); row<_n_rows; row++) {
		for (std::size_t col(0); col<_n_cols; col++) {
			for (std::size_t new_row(row*row_blk_size); new_row<(row+1)*row_blk_size; new_row++) {
				for (std::size_t new_col(col*row_blk_size); new_col<(col+1)*col_blk_size; new_col++) {
					new_raster_data[new_row*n_cols+new_col] = _raster_data[row*_n_cols+col];
				}
			}
		}
	}

	std::swap(_raster_data, new_raster_data);
	_cell_size /= row_blk_size;
	_n_cols = n_cols;
	_n_rows = n_rows;

	delete [] new_raster_data;
}

Raster::~Raster()
{
	if (_raster_data != NULL)
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

void Raster::writeRasterAsASC(std::ostream &os) const
{
	// write header
	os << "ncols " << _n_rows << std::endl;
	os << "nrows " << _n_cols << std::endl;
	os << "xllcorner " << _ll_pnt[0] << std::endl;
	os << "yllcorner " << _ll_pnt[1] << std::endl;
	os << "cellsize " <<  _cell_size << std::endl;
	os << "NODATA_value " << _no_data_val << std::endl;

	// write data
	for (unsigned row(0); row<_n_rows; row++) {
		for (unsigned col(0); col<_n_cols; col++) {
			os << _raster_data[(_n_rows-row-1)*_n_cols+col] << " ";
		}
		os << std::endl;
	}
}

Raster* Raster::getRasterFromSurface(Surface const& sfc, double cell_size, double no_data_val)
{
	Point const& ll (sfc.getAABB().getMinPoint());
	Point const& ur (sfc.getAABB().getMaxPoint());

	std::size_t n_cols = static_cast<size_t>(fabs(ur[0]-ll[0]) / cell_size)+1;
	std::size_t n_rows = static_cast<size_t>(fabs(ur[1]-ll[1]) / cell_size)+1;
	const size_t n_triangles (sfc.getNTriangles());
	double *z_vals (new double[n_cols*n_rows]);
	if (!z_vals) {
		std::cout << "DEBUG: CreateRaster::getRaster " << n_cols << " x " << n_rows << " to big" << std::endl;
		return NULL;
	}
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

	return new Raster(n_cols, n_rows, ll[0], ll[1], cell_size, -9999, z_vals);
}

Raster* Raster::getRasterFromASCFile(std::string const& fname)
{
	std::ifstream in(fname.c_str());

	if (!in.is_open()) {
		std::cout << "Raster::getRasterFromASCFile() - Could not open file..." << fname << std::endl;
		return NULL;
	}

	// header information
	std::size_t n_cols(0), n_rows(0);
	double xllcorner(0.0), yllcorner(0.0), cell_size(0.0), no_data_val(-9999);

	if (readASCHeader(in, n_cols, n_rows, xllcorner, yllcorner, cell_size, no_data_val)) {
		double* values = new double[n_cols*n_rows];
		std::string s;
		// read the data into the double-array
		for (size_t j(0); j < n_rows; ++j) {
			size_t idx ((n_rows - j - 1) * n_cols);
			for (size_t i(0); i < n_cols; ++i) {
				in >> s;
				values[idx+i] = strtod(BaseLib::replaceString(",", ".", s).c_str(),0);

			}
		}
		in.close();
		return new Raster(n_cols, n_rows, xllcorner, yllcorner,
						cell_size, no_data_val, values);
	} else {
		std::cout << "Raster::getRasterFromASCFile() - could not read header of file " << fname << std::endl;
		return NULL;
	}
}

bool Raster::readASCHeader(std::ifstream &in, size_t &n_cols, std::size_t &n_rows,
				double &xllcorner, double &yllcorner, double &cell_size, double &no_data_val)
{
	std::string tag, value;

	in >> tag;
	if (tag.compare("ncols") == 0) {
		in >> value;
		n_cols = atoi(value.c_str());
	} else return false;

	in >> tag;
	if (tag.compare("nrows") == 0) {
		in >> value;
		n_rows = atoi(value.c_str());
	} else return false;

	in >> tag;
	if (tag.compare("xllcorner") == 0) {
		in >> value;
		xllcorner = strtod(BaseLib::replaceString(",", ".", value).c_str(), 0);
	} else return false;

	in >> tag;
	if (tag.compare("yllcorner") == 0) {
		in >> value;
		yllcorner = strtod(BaseLib::replaceString(",", ".", value).c_str(), 0);
	} else return false;

	in >> tag;
	if (tag.compare("cellsize") == 0) {
		in >> value;
		cell_size = strtod(BaseLib::replaceString(",", ".", value).c_str(), 0);
	} else return false;

	in >> tag;
	if (tag.compare("NODATA_value") == 0) {
		in >> value;
		no_data_val = strtod(BaseLib::replaceString(",", ".", value).c_str(), 0);
	} else return false;

	return true;
}

} // end namespace GeoLib
