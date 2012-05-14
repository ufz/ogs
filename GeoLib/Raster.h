/*
 * Raster.h
 *
 *  Created on: Sep 7, 2011
 *      Author: TF
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
	double* getRasterFromSurface (Surface const& sfc, size_t &n_x_pnts, size_t &n_y_pnts) const;
	virtual ~Raster();
private:
	double _cell_size;
	double _no_data_val;
};

}

#endif /* RASTER_H_ */
