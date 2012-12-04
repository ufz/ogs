/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkRaster.h
 *
 * Created on 2012-02-01 by Karsten Rink
 */
#ifndef VTKRASTER_H
#define VTKRASTER_H

#include <string>

class vtkImageAlgorithm;
class vtkImageImport;
class vtkImageReader2;

/**
 * \brief Loading raster data such as images or ArcGIS-data into VTK image data structures.
 *
 * The VtkRaster class enables loading of raster data such as images or ArcGIS-data. Supported image formats are .....
 * Georeferenced data can be imported via the GeoTIFF- or asc-format.
 */
class VtkRaster
{

public:
	/**
	 * \brief Loads an image- or raster-file into an vtkImageAlgorithm-Object.
	 *
	 * Public method for loading all data formats. Internally the method automatically differentiates between
	 * images and georeferenced files and then calls the appropriate method for reading the file.
	 * \param fileName Filename of the file that should be loaded.
	 * \param x0 X-coordinate of the upper left corner of the data set, the default value is 0.
	 * \param y0 Y-coordinate of the upper left corner of the data set, the default value is 0.
	 * \param delta The size of each pixel in the image which is needed for correctly displaying the data, the default value is 1.
	 * \return The ImageAlgorithm-object.
	 */
    static vtkImageAlgorithm* loadImage(const std::string &fileName,
                                        double& x0,
                                        double& y0,
                                        double& delta);

	/**
	 * \brief Returns a VtkImageAlgorithm from an array of pixel values and some image meta data.
	 */
    static vtkImageImport* loadImageFromArray(double* data_array,
											  double &x0,
											  double &y0,
											  unsigned &width,
											  unsigned &height,
											  double &delta,
											  double no_data = -9999);

private:
	/**
	 * Loads ArcGIS asc-files to a QPixmap object and automatically does a contrast stretching to adjust values to 8 bit greyscale images.
	 * \param fileName Filename of the file that should be loaded.
	 * \param raster The QPixmap into which the raster data will be written.
	 * \param origin The upper left corner of the data set
	 * \param delta The size of each pixel in the image which is needed for correctly displaying the data.
	 * \return A vtkImageImport-object (derived from vtkImageAlgorithm).
	 */
#ifdef libgeotiff_FOUND
	static vtkImageImport* loadImageFromTIFF(const std::string &fileName,
											 double& x0, double& y0, double& delta);
#endif

	/**
	 * Loads image files into a QPixmap object. Since images are not geo-referenced no origin point will be returned.
	 * \param fileName Filename of the file that should be loaded.
	 * \return vtkImageReader2-object containing the image data.
	 */
	static vtkImageReader2* loadImageFromFile(const std::string &fileName);

	/// Converts an uint32-number into a quadruple representing RGBA-colours for a pixel.
	static void uint32toRGBA(const unsigned int s, int* p);
};

#endif //VTKRASTER_H
