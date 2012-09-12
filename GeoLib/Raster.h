/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
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
	Raster(double cell_size=1, double no_data_val=9999);
	void setCellSize(double cell_size);
	void setNoDataVal (double no_data_val);
	double* getRasterFromSurface (Surface const& sfc, std::size_t &n_x_pnts, std::size_t &n_y_pnts) const;
	virtual ~Raster();
private:
	double _cell_size;
	double _no_data_val;
};

}

#endif /* RASTER_H_ */
