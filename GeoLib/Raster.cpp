/**
 * @file Raster.cpp
 * @author Thomas Fischer
 * @date 2011-09-07
 * @brief Implementation of the GeoLib::Raster class.
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

double Raster::getValueAtPoint(const GeoLib::Point &pnt)
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

void Raster::writeRasterAsASC(std::ostream &os) const
{
	// write header
	os << "ncols " << _n_cols << "\n";
	os << "nrows " << _n_rows << "\n";
	os << "xllcorner " << _ll_pnt[0] << "\n";
	os << "yllcorner " << _ll_pnt[1] << "\n";
	os << "cellsize " <<  _cell_size << "\n";
	os << "NODATA_value " << _no_data_val << "\n";

	// write data
	for (unsigned row(0); row<_n_rows; row++) {
		for (unsigned col(0); col<_n_cols; col++) {
			os << _raster_data[(_n_rows-row-1)*_n_cols+col] << " ";
		}
		os << "\n";
	}
}

Raster* Raster::readRaster(std::string const& fname)
{
	std::string ext (BaseLib::getFileExtension(fname));
	std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
	if (ext.compare("asc") == 0)
		return getRasterFromASCFile(fname);
	if (ext.compare("grd") == 0)
		return getRasterFromSurferFile(fname);
	return nullptr;
}

Raster* Raster::getRasterFromASCFile(std::string const& fname)
{
	std::ifstream in(fname.c_str());

	if (!in.is_open()) {
		WARN("Raster::getRasterFromASCFile(): Could not open file %s.", fname.c_str());
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
			const size_t idx ((n_rows - j - 1) * n_cols);
			for (size_t i(0); i < n_cols; ++i) {
				in >> s;
				values[idx+i] = strtod(BaseLib::replaceString(",", ".", s).c_str(),0);

			}
		}
		in.close();
		Raster *raster(new Raster(n_cols, n_rows, xllcorner, yllcorner,
						cell_size, values, values+n_cols*n_rows, no_data_val));
		delete [] values;
		return raster;
	} else {
		WARN("Raster::getRasterFromASCFile(): Could not read header of file %s", fname.c_str());
		return NULL;
	}
}

bool Raster::readASCHeader(std::ifstream &in, std::size_t &n_cols, std::size_t &n_rows,
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

Raster* Raster::getRasterFromSurferFile(std::string const& fname)
{
	std::ifstream in(fname.c_str());

	if (!in.is_open()) {
		ERR("Raster::getRasterFromSurferFile() - Could not open file %s", fname.c_str());
		return NULL;
	}

	// header information
	std::size_t n_cols(0), n_rows(0);
	double xllcorner(0.0), yllcorner(0.0), cell_size(0.0), min(0.0), max(0.0);

	if (readSurferHeader(in, n_cols, n_rows, xllcorner, yllcorner, cell_size, min, max)) 
	{
		const double no_data_val (min-1);
		double* values = new double[n_cols*n_rows];
		std::string s;
		// read the data into the double-array
		for (size_t j(0); j < n_rows; ++j) 
		{
			const size_t idx (j * n_cols);
			for (size_t i(0); i < n_cols; ++i) 
			{
				in >> s;
				const double val (strtod(BaseLib::replaceString(",", ".", s).c_str(),0));
				values[idx+i] = (val > max || val < min) ? no_data_val : val;
			}
		}
		in.close();
		Raster *raster(new Raster(n_cols, n_rows, xllcorner, yllcorner,
						cell_size, values, values+n_cols*n_rows, no_data_val));
		delete [] values;
		return raster;
	} else {
		ERR("Raster::getRasterFromASCFile() - could not read header of file %s", fname.c_str());
		return NULL;
	}
}

bool Raster::readSurferHeader(std::ifstream &in, size_t &n_cols, std::size_t &n_rows,
				double &xllcorner, double &yllcorner, double &cell_size, double &min, double &max)
{
	std::string tag;
	
	in >> tag;

	if (tag.compare("DSAA") != 0)
	{
		ERR("Error in readSurferHeader() - No Surfer file.");
		return false;
	}
	else
	{
		in >> n_cols >> n_rows;
		in >> min >> max;
		xllcorner = min;
		cell_size = (max-min)/static_cast<double>(n_cols);

		in >> min >> max;
		yllcorner = min;

		if (ceil((max-min)/static_cast<double>(n_rows)) == ceil(cell_size))
			cell_size = ceil(cell_size);
		else
		{
			ERR("Error in readSurferHeader() - Anisotropic cellsize detected.");
			return 0;
		}
		in >> min >> max;
	}

	return true;
}

} // end namespace GeoLib
