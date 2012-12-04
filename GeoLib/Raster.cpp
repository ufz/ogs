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
				for (std::size_t new_col(col*col_blk_size); new_col<(col+1)*col_blk_size; new_col++) {
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

Raster* Raster::getRasterFromSurface(GeoLib::Surface const& sfc, double cell_size, double no_data_val)
{
	GeoLib::Point const& ll (sfc.getAABB().getMinPoint());
	GeoLib::Point const& ur (sfc.getAABB().getMaxPoint());

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

	return new Raster(n_cols, n_rows, ll[0], ll[1], cell_size, -9999, z_vals);
}

Raster* Raster::getRasterFromASCFile(std::string const& fname)
{
	unsigned width(0), height(0);
	double x0(0), y0(0), delta(1);
	double no_data(-9999);

	double* data = loadDataFromASC(fname, x0, y0, height, width, delta, no_data);

	if (data)
		return new Raster(width, height, x0, y0, delta, no_data, data);
	else
		return nullptr;
}

bool Raster::readASCHeader(ascHeader &header, std::ifstream &in)
{
	std::string line, tag, value;

	in >> tag;
	if (tag.compare("ncols") == 0)
	{
		in >> value;
		header.ncols = atoi(value.c_str());
	}
	else
		return false;
	in >> tag;
	if (tag.compare("nrows") == 0)
	{
		in >> value;
		header.nrows = atoi(value.c_str());
	}
	else
		return false;
	in >> tag;
	if (tag.compare("xllcorner") == 0)
	{
		in >> value;
		header.x = strtod(BaseLib::replaceString(",", ".", value).c_str(),0);
	}
	else
		return false;
	in >> tag;
	if (tag.compare("yllcorner") == 0)
	{
		in >> value;
		header.y = strtod(BaseLib::replaceString(",", ".", value).c_str(),0);
	}
	else
		return false;
	in >> tag;
	if (tag.compare("cellsize") == 0)
	{
		in >> value;
		header.cellsize = strtod(BaseLib::replaceString(",", ".", value).c_str(),0);
	}
	else
		return false;
	in >> tag;
	if (tag.compare("NODATA_value") == 0)
	{
		in >> value;
		header.noData = value.c_str();
	}
	else
		return false;

	// correct raster position by half a pixel for correct visualisation
	// argh! wrong! correction has to happen in visualisation object, otherwise the actual data is wrong
	//header.x = header.x + (header.cellsize / 2);
	//header.y = header.y + (header.cellsize / 2);

	return true;
}

double* Raster::loadDataFromASC(const std::string &fileName,
                                   double &x0,
                                   double &y0,
                                   unsigned &width,
                                   unsigned &height,
                                   double &delta,
								   double &no_data)
{
	std::ifstream in( fileName.c_str() );

	if (!in.is_open())
	{
		std::cout << "VtkRaster::loadImageFromASC() - Could not open file..." << std::endl;
		return NULL;
	}

	ascHeader header;

	if (readASCHeader(header, in))
	{
		x0     = header.x;
		y0     = header.y;
		width  = header.ncols;
		height = header.nrows;
		delta  = header.cellsize;

		double* values = new double[header.ncols * header.nrows];

		int col_index(0);
		int noData = atoi(header.noData.c_str());
		std::string s("");
		// read the file into a double-array
		for (int j = 0; j < header.nrows; ++j)
		{
			col_index = (header.nrows - j - 1) * header.ncols;
			for (int i = 0; i < header.ncols; ++i)
			{
				in >> s;
				unsigned index = col_index+i;
				values[index] = strtod(BaseLib::replaceString(",", ".", s).c_str(),0);
			}
		}

		in.close();
		return values;
	}
	return nullptr;
}

bool Raster::readSurferHeader(ascHeader &header, std::ifstream &in)
{
	std::string line, tag, value;
	double min, max;

	in >> tag;

	if (tag.compare("DSAA") != 0)
	{
		std::cout << "Error in readSurferHeader() - No Surfer file..." << std::endl;
		return false;
	}
	else
	{
		in >> header.ncols >> header.nrows;
		in >> min >> max;
		header.x = min;
		header.cellsize = (max-min)/(double)header.ncols;

		in >> min >> max;
		header.y = min;

		if (ceil((max-min)/(double)header.nrows) == ceil(header.cellsize))
			header.cellsize = ceil(header.cellsize);
		else
		{
			std::cout << "Error in readSurferHeader() - Anisotropic cellsize detected..." << std::endl;
			return 0;
		}
		in >> min >> max; // ignore min- and max-values

		header.noData = "1.70141E+038";
	}

	return true;
}

double* Raster::loadDataFromSurfer(const std::string &fileName,
                                   double &x0,
                                   double &y0,
                                   unsigned &width,
                                   unsigned &height,
                                   double &delta,
								   double &no_data)
{
	std::ifstream in( fileName.c_str() );

	if (!in.is_open())
	{
		std::cout << "VtkRaster::loadImageFromSurfer() - Could not open file..." << std::endl;
		return NULL;
	}

	ascHeader header;

	if (readSurferHeader(header, in))
	{
		x0     = header.x;
		y0     = header.y;
		width  = header.ncols;
		height = header.nrows;
		delta  = header.cellsize;

		double* values = new double[header.ncols * header.nrows];

		int col_index(0);
		int noData = -9999;
		std::string s("");
		// read the file into a double-array
		for (int j = 0; j < header.nrows; ++j)
		{
			col_index = j * header.ncols;
			for (int i = 0; i < header.ncols; ++i)
			{
				in >> s;
				if (s.compare(header.noData) == 0)
					s = "-9999";
				unsigned index = col_index+i;
				values[index] = strtod(BaseLib::replaceString(",", ".", s).c_str(),0);
			}
		}

		in.close();
		return values;
	}
	return nullptr;
}

