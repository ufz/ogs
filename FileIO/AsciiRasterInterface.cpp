/**
 * @file AsciiRasterInterface.cpp
 * @author Karsten Rink
 * @date 2014-09-10
 * @brief Implementation of the AsciiRasterInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "AsciiRasterInterface.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "Raster.h"

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

namespace FileIO {

GeoLib::Raster* AsciiRasterInterface::readRaster(std::string const& fname)
{
	std::string ext (BaseLib::getFileExtension(fname));
	std::transform(ext.begin(), ext.end(), ext.begin(), tolower);
	if (ext.compare("asc") == 0)
		return getRasterFromASCFile(fname);
	if (ext.compare("grd") == 0)
		return getRasterFromSurferFile(fname);
	return nullptr;
}

GeoLib::Raster* AsciiRasterInterface::getRasterFromASCFile(std::string const& fname)
{
	std::ifstream in(fname.c_str());

	if (!in.is_open()) {
		WARN("Raster::getRasterFromASCFile(): Could not open file %s.", fname.c_str());
		return nullptr;
	}

	// header information
	std::size_t n_cols(0), n_rows(0);
	double xllcorner(0.0), yllcorner(0.0), cell_size(0.0), no_data_val(-9999);

	if (readASCHeader(in, n_cols, n_rows, xllcorner, yllcorner, cell_size, no_data_val)) {
		double* values = new double[n_cols*n_rows];
		std::string s;
		// read the data into the double-array
		for (std::size_t j(0); j < n_rows; ++j) {
			const size_t idx ((n_rows - j - 1) * n_cols);
			for (size_t i(0); i < n_cols; ++i) {
				in >> s;
				values[idx+i] = strtod(BaseLib::replaceString(",", ".", s).c_str(),0);

			}
		}
		in.close();
		GeoLib::Raster *raster(new GeoLib::Raster(n_cols, n_rows, xllcorner, yllcorner,
						cell_size, values, values+n_cols*n_rows, no_data_val));
		delete [] values;
		return raster;
	} else {
		WARN("Raster::getRasterFromASCFile(): Could not read header of file %s", fname.c_str());
		return nullptr;
	}
}

bool AsciiRasterInterface::readASCHeader(std::ifstream &in, std::size_t &n_cols, std::size_t &n_rows,
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

GeoLib::Raster* AsciiRasterInterface::getRasterFromSurferFile(std::string const& fname)
{
	std::ifstream in(fname.c_str());

	if (!in.is_open()) {
		ERR("Raster::getRasterFromSurferFile() - Could not open file %s", fname.c_str());
		return nullptr;
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
		GeoLib::Raster *raster(new GeoLib::Raster(n_cols, n_rows, xllcorner, yllcorner, cell_size,
						                          values, values+n_cols*n_rows, no_data_val));
		delete [] values;
		return raster;
	} else {
		ERR("Raster::getRasterFromASCFile() - could not read header of file %s", fname.c_str());
		return nullptr;
	}
}

bool AsciiRasterInterface::readSurferHeader(std::ifstream &in, size_t &n_cols, std::size_t &n_rows,
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

void AsciiRasterInterface::writeRasterAsASC(GeoLib::Raster const& raster, std::string const& file_name)
{
    GeoLib::Point const& origin (raster.getOrigin());
    unsigned const nCols (raster.getNCols());
    unsigned const nRows (raster.getNRows());

    // write header
    std::ofstream out(file_name);
    out << "ncols " << nCols << "\n";
    out << "nrows " << nRows << "\n";
    out << "xllcorner " << origin[0] << "\n";
    out << "yllcorner " << origin[1] << "\n";
    out << "cellsize " <<  raster.getRasterPixelSize() << "\n";
    out << "NODATA_value " << raster.getNoDataValue() << "\n";

    // write data
    double const*const elevation(raster.begin());
    for (unsigned row(0); row < nRows; ++row) 
    {
        for (unsigned col(0); col < nCols; ++col)
        {
            out << elevation[(nRows-row-1) * nCols + col] << " ";
        }
        out << "\n";
    }
    out.close();
}

} // end namespace GeoLib
