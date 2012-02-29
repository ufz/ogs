/**
 * \file VtkRaster.h
 * 2012/02/01 KR Initial implementation
 *
 * (Re-implementation of the old OGSRaster class.)
 */
#ifndef VTKRASTER_H
#define VTKRASTER_H

#include <fstream>

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
	 * \brief Loads an ASC file into a double array
	 *
	 * \param fileName Filename of the file that should be loaded.
	 * \param x0 The x-coordinate of the origin.
	 * \param y0 The y-coordinate of the origin.
	 * \param width The width of the image.
	 * \param height The height of the image
	 * \param delta The size of each pixel in the image which is needed for correctly displaying the data.
	 * \return A float-array of pixel values incl. opacity (noData values are transparent)
	 */
	static float* loadDataFromASC(const std::string &fileName,
	                              double &x0,
	                              double &y0,
	                              size_t &width,
	                              size_t &height,
	                              double &delta);

	static float* loadDataFromSurfer(const std::string &fileName,
	                              double &x0,
	                              double &y0,
	                              size_t &width,
	                              size_t &height,
	                              double &delta);

    static vtkImageImport* loadImageFromArray(double* data_array, 
											  double &x0, 
											  double &y0, 
											  size_t &width,
											  size_t &height,
											  double &delta,
											  double noData);

private:
	/**
	 * Loads ArcGIS asc-files to a vtkImageImport object.
	 * \param fileName Filename of the file that should be loaded.
	 * \param x0 The x-coordinate of the origin.
	 * \param y0 The y-coordinate of the origin.
	 * \param delta The size of each pixel in the image which is needed for correctly displaying the data.
	 * \return A vtkImageImport-object (derived from vtkImageAlgorithm).
	 */
    static vtkImageImport* loadImageFromASC(const std::string &fileName,
                                            double& x0, double& y0, double& delta);

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

	/// Converts an uint32-number into a quadruple representing RGBA-colours for a pixel.
	static void uint32toRGBA(const unsigned int s, int* p);
};

#endif //VTKRASTER_H
