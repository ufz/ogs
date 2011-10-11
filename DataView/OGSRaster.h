/**
 * \file OGSRaster.h
 * 11/01/2010 KR Initial implementation
 *
 */
#ifndef OGSRASTER_H
#define OGSRASTER_H

#include <fstream>

class QImage;
class QPointF;
class QString;

/**
 * \brief Loading of raster data such as images or ArcGIS-data.
 *
 * The OGSRaster class enables loading of raster data such as images or ArcGIS-data. Supported image formats are specified
 * by the Qt class QPixmap. Georeferenced data can be imported via the GeoTIFF- or asc-format .
 */
class OGSRaster
{
	/// Data structure for the asc-file header.
	struct ascHeader
	{
		int ncols;
		int nrows;
		double x;
		double y;
		double cellsize;
		int noData;
	};

public:
	/**
	 * \brief Loads an image- or raster-file into an QImage.
	 *
	 * Public method for loading all data formats. Internally the method automatically differentiates between
	 * images and georeferenced files and then calls the appropriate method for reading the file.
	 * \param fileName Filename of the file that should be loaded.
	 * \param raster The QImage into which the raster data will be written.
	 * \param origin The upper left corner of the data set, the default value is (0,0).
	 * \param scalingFactor The size of each pixel in the image which is needed for re-scaling the data, the default value is 1.
	 * \param autoscale Determines if the histogram of the raster will be contrast-stretched to [0, 255]. If false, the streching process will be skipped.
	 * \param mirrorX Mirror around x-axis.
	 * \return True if the raster data was loaded correctly, false otherwise.
	 */
	static bool loadImage(const QString &fileName,
	                      QImage &raster,
	                      QPointF &origin,
	                      double &scalingFactor,
	                      bool autoscale = true,
	                      bool mirrorX = false);

	/**
	 * \brief Loads an ASC file into a double array
	 *
	 * \param fileName Filename of the file that should be loaded.
	 * \param x0 The x-coordinate of the origin.
	 * \param y0 The y-coordinate of the origin.
	 * \param width The width of the image.
	 * \param height The height of the image
	 * \param delta The size of each pixel in the image which is needed for re-scaling the data.
	 * \return True if the raster data was loaded correctly, false otherwise.
	 */
	static double* loadDataFromASC(const QString &fileName,
	                               double& x0,
	                               double &y0,
	                               size_t &width,
	                               size_t &height,
	                               double &delta);

	/// Converts raster to an 8 bit greyscale image that is contrast-stretched in [min:max].
	static void convertToGreyscale(QImage &raster, const int &min = 0, const int &max = 255);

	/// Returns an int-array containing the raster converted an 8 bit greyscale values that are contrast-stretched in [min:max].
	static int* getGreyscaleData(QImage &raster, const int &min = 0, const int &max = 255);

	/// Returns the maximum intensity of raster.
	static int getMaxValue(const QImage &raster);

	/// Returns the minimum intensity of raster.
	static int getMinValue(const QImage &raster);

private:
	/**
	 * Loads ArcGIS asc-files to a QPixmap object and automatically does a contrast stretching to adjust values to 8 bit greyscale images.
	 * \param fileName Filename of the file that should be loaded.
	 * \param raster The QPixmap into which the raster data will be written.
	 * \param origin The upper left corner of the data set
	 * \param scalingFactor
	 * \param autoscale
	 * \return True if the raster data was loaded correctly, false otherwise.
	 */
	static bool loadImageFromASC(const QString &fileName,
	                             QImage &raster,
	                             QPointF &origin,
	                             double &scalingFactor,
	                             bool autoscale = true);

	/**
	 * Loads ArcGIS asc-files to a QPixmap object and automatically does a contrast stretching to adjust values to 8 bit greyscale images.
	 * \param fileName Filename of the file that should be loaded.
	 * \param raster The QPixmap into which the raster data will be written.
	 * \param origin The upper left corner of the data set
	 * \param scalingFactor
	 * \return True if the raster data was loaded correctly, false otherwise.
	 */
#ifdef libgeotiff_FOUND
	static bool loadImageFromTIFF(const QString &fileName,
	                              QImage &raster,
	                              QPointF &origin,
	                              double &scalingFactor);
#endif

	/**
	 * Loads image files into a QPixmap object. Since images are not geo-referenced no origin point will be returned.
	 * \param fileName Filename of the file that should be loaded.
	 * \param raster The QPixmap into which the raster data will be written.
	 * \return True if the raster data was loaded correctly, false otherwise.
	 */
	static bool loadImageFromFile(const QString &fileName, QImage &raster);

	/**
	 * Reads the header of an ArcGIS asc-file.
	 * \param header The ascHeader-object into which all the information will be written.
	 * \param in FileInputStream used for reading the data.
	 * \return True if the header could be read correctly, false otherwise.
	 */
	static bool readASCHeader(ascHeader &header, std::ifstream &in);
	static void uint32toRGBA(const unsigned int s, int* p);
};

#endif //OGSRASTER_H
