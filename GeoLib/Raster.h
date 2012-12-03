/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Raster.h
 *
 * Created on 2011-09-07 by Thomas Fischer
 */

#ifndef RASTER_H_
#define RASTER_H_

#include "Surface.h"

namespace GeoLib {

class Raster {
public:
	Raster(std::size_t n_cols, std::size_t n_rows, double xllcorner, double yllcorner,
					double cell_size = 1, double no_data_val = -9999, double* raster_data = NULL);

	std::size_t getNCols() const { return _n_cols; }
	std::size_t getNRows() const { return _n_rows; }

	/**
	 * get the distance between raster pixels
	 * @return
	 */
	double getRasterPixelDistance() const { return _cell_size; }

	/**
	 * get the origin of lower left raster cell
	 * @return the origin of the raster
	 */
	GeoLib::Point const& getOrigin() const;

	void refineRaster(std::size_t n_cols, std::size_t n_rows);
	double const* getRasterData() const { return _raster_data; }
	virtual ~Raster();

	void writeRasterAsASC(std::ostream &os) const;

	static Raster* getRasterFromASCFile(std::string const& fname);
	static Raster* getRasterFromSurface(Surface const& sfc, double cell_size, double no_data_val = -9999);

private:
	static bool readASCHeader(std::ifstream &in, std::size_t &n_cols, std::size_t &n_rows,
					double &xllcorner, double &yllcorner, double &cell_size, double &no_data_val);
	void setCellSize(double cell_size);
	void setNoDataVal (double no_data_val);

	std::size_t _n_cols;
	std::size_t _n_rows;
	GeoLib::Point _ll_pnt;
	double _cell_size;
	double _no_data_val;
	double* _raster_data;
};

}

#endif /* RASTER_H_ */
