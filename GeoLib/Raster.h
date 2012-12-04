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

	/**
	 * \brief Loads an ASC file into a double array.
	 * The array alternates between pixel values and their respective alpha-values, i.e.
	 * result = { pixel0-value; pixel0-alpha, pixel1-value; pixel1-alpha; ... }
	 *
	 * \param fileName Filename of the file that should be loaded.
	 * \param x0 The x-coordinate of the origin.
	 * \param y0 The y-coordinate of the origin.
	 * \param width The width of the image.
	 * \param height The height of the image
	 * \param delta The size of each pixel in the image which is needed for correctly displaying the data.
	 * \param no_data The value that signifies that no meaningful data is given at certain pixels
	 * \return A float-array of pixel values incl. opacity (noData values are transparent)
	 */
	static double* loadDataFromASC(const std::string &fileName,
	                              double &x0,
	                              double &y0,
	                              unsigned &width,
	                              unsigned &height,
	                              double &delta,
								  double &no_data);

	/**
	 * \brief Loads a Surfer file into a double array.
	 * Works exactly like loadDataFromASC().
	 */
	static double* loadDataFromSurfer(const std::string &fileName,
	                              double &x0,
	                              double &y0,
	                              unsigned &width,
	                              unsigned &height,
	                              double &delta,
								  double &no_data);

	static Raster* getRasterFromASCFile(std::string const& fname);
	static Raster* getRasterFromSurface(GeoLib::Surface const& sfc, double cell_size, double no_data_val = -9999);

private:
	void setCellSize(double cell_size);
	void setNoDataVal (double no_data_val);

	/// Data structure for the asc-file header.
	struct ascHeader
	{
		int ncols;
		int nrows;
		double x;
		double y;
		double cellsize;
		std::string noData;
	};

	/**
	 * Reads the header of an ArcGIS asc-file.
	 * \param header The ascHeader-object into which all the information will be written.
	 * \param in FileInputStream used for reading the data.
	 * \return True if the header could be read correctly, false otherwise.
	 */
	static bool readASCHeader(ascHeader &header, std::ifstream &in);

	/**
	 * Reads the header of a Surfer grd-file.
	 * \param header The ascHeader-object into which all the information will be written.
	 * \param in FileInputStream used for reading the data.
	 * \return True if the header could be read correctly, false otherwise.
	 */
	static bool readSurferHeader(ascHeader &header, std::ifstream &in);

	std::size_t _n_cols;
	std::size_t _n_rows;
	GeoLib::Point _ll_pnt;
	double _cell_size;
	double _no_data_val;
	double* _raster_data;
};


#endif /* RASTER_H_ */
