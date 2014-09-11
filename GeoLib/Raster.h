/**
 * @file Raster.h
 * @author Thomas Fischer
 * @date 2011-09-07
 * @brief Definition of the GeoLib::Raster class.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef RASTER_H_
#define RASTER_H_

#include "Surface.h"

namespace GeoLib {

/**
 * @brief Class Raster is used for managing raster data.
 *
 * A raster consists of the meta data like number of columns and rows, the lower
 * left point, the size of a raster pixel and a value for invalid data pixels.
 * Additional the object needs the raster data itself. The raster data will be
 * copied from the constructor. The destructor will release the memory.
 */
class Raster {
public:
	typedef double const* const_iterator;
	typedef double* iterator;

	/**
	 * @brief Constructor for an object of class Raster. The raster data have
	 * to be handed over via input iterators. Deploying iterators has the
	 * advantage that the use of the class is independent from the input
	 * container.
	 * @param n_cols number of columns
	 * @param n_rows number of rows
	 * @param xllcorner the \f$x\f$ coordinate of lower left point
	 * @param yllcorner the \f$y\f$ coordinate of lower left point
	 * @param cell_size the size of a raster pixel
	 * @param begin input iterator pointing to the first element of the data
	 * @param end input iterator pointing to the last element of the data, end have to be reachable from begin
	 * @param no_data_val value for raster pixels that have not valid values
	 */
	template<typename InputIterator>
	Raster(std::size_t n_cols, std::size_t n_rows, double xllcorner, double yllcorner,
					double cell_size, InputIterator begin, InputIterator end, double no_data_val = -9999) :
		_n_cols(n_cols), _n_rows(n_rows), _ll_pnt(xllcorner, yllcorner, 0.0), _cell_size(cell_size),
		_no_data_val(no_data_val), _raster_data(new double[n_cols*n_rows])
	{
		iterator raster_it(_raster_data);
		for (InputIterator it(begin); it != end; ++it) {
			*raster_it = *it;
			raster_it++;
		}
	}

	/**
	 * get the number of columns for the raster
	 */
	std::size_t getNCols() const { return _n_cols; }
	/**
	 * get the number of rows for the raster
	 */
	std::size_t getNRows() const { return _n_rows; }

	/**
	 * get the distance between raster pixels
	 */
	double getRasterPixelSize() const { return _cell_size; }

	/**
	 * get the origin of lower left raster cell
	 * @return the origin of the raster
	 */
	GeoLib::Point const& getOrigin() const;

	/**
	 * Refine the raster using scaling as a refinement parameter.
	 */
	void refineRaster(std::size_t scaling);

	double getNoDataValue() const { return _no_data_val; }

	/**
	 * Constant iterator that is pointing to the first raster pixel value.
	 * @return constant iterator
	 */
	const_iterator begin() const { return _raster_data; }
	/**
	 * Constant iterator that is pointing to the last raster pixel value.
	 * @return constant iterator
	 */
	const_iterator end() const { return _raster_data + _n_rows*_n_cols; }

	/**
	 * Returns the raster value at the position of the given point.
	 */
	double getValueAtPoint(const GeoLib::Point &pnt) const;

	/// Interpolates the elevation of the given point based on the 8-neighbourhood of the raster cell it is located on
	double interpolateValueAtPoint(const GeoLib::Point &pnt) const;
    
	/// Checks if the given point is located within the (x,y)-extension of the raster.
	bool isPntOnRaster(GeoLib::Point const& node) const;

	~Raster();

	/// Creates a Raster based on a GeoLib::Surface
	static Raster* getRasterFromSurface(Surface const& sfc, double cell_size, double no_data_val = -9999);

private:
	void setCellSize(double cell_size);
	void setNoDataVal (double no_data_val);

	/**
	 * number of columns of the raster
	 */
	std::size_t _n_cols;
	/**
	 * number of rows of the raster
	 */
	std::size_t _n_rows;
	/**
	 * lower left point of the raster
	 */
	GeoLib::Point _ll_pnt;
	/**
	 * cell size - each cell is quadratic
	 */
	double _cell_size;
	/**
	 * value for cells there is no data for
	 */
	double _no_data_val;
	/**
	 * raw raster data
	 */
	double* _raster_data;
};

}

#endif /* RASTER_H_ */
