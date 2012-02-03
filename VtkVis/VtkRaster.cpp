/**
 * \file VtkRaster.cpp
 * 2012/02/01 KR Initial implementation
 */

#include "VtkRaster.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <QFileInfo>

#include "OGSError.h"
#include "StringTools.h"

#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vtkImageAlgorithm.h>
#include <vtkImageImport.h>


/*
#ifdef libgeotiff_FOUND
#include "geo_tiffp.h"
#include "xtiffio.h"
#endif
*/
vtkImageAlgorithm* VtkRaster::loadImage(const std::string &fileName,
                                        double& x0, double& y0,
                                        double& delta, bool autoscale)
{
	QFileInfo fileInfo(QString::fromStdString(fileName));

	if (fileInfo.suffix().toLower() == "asc")
	{
        return loadImageFromASC(fileName, x0, y0, delta, autoscale);
	}
/*
#ifdef libgeotiff_FOUND
	else if (fileInfo.suffix().toLower() == "tif")
	{
		if (!loadImageFromTIFF(fileName, raster))
			return false;
	}
#endif
	else if (!loadImageFromFile(fileName, raster))
		return false;
*/
}

vtkImageImport* VtkRaster::loadImageFromASC(const std::string &fileName,
                                            double& x0, double& y0,
                                            double& delta, bool autoscale)
{
	size_t width(0), height(0);
	float* data = loadDataFromASC(fileName, x0, y0, width, height, delta);

	vtkImageImport* image = vtkImageImport::New();
		image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
		image->SetDataExtent(0, width-1, 0, height-1, 0, 0);
		image->SetDataExtentToWholeExtent();
		image->SetDataScalarTypeToFloat();
		image->SetNumberOfScalarComponents(1);
		image->SetImportVoidPointer(data, 0);
		image->Update();

	return image;
}

bool VtkRaster::readASCHeader(ascHeader &header, std::ifstream &in)
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
		header.x = strtod(replaceString(",", ".", value).c_str(),0);
	}
	else
		return false;
	in >> tag;
	if (tag.compare("yllcorner") == 0)
	{
		in >> value;
		header.y = strtod(replaceString(",", ".", value).c_str(),0);
	}
	else
		return false;
	in >> tag;
	if (tag.compare("cellsize") == 0)
	{
		in >> value;
		header.cellsize = strtod(replaceString(",", ".", value).c_str(),0);
	}
	else
		return false;
	in >> tag;
	if (tag.compare("NODATA_value") == 0)
	{
		in >> value;
		header.noData = atoi(value.c_str());
	}
	else
		return false;

	// correct raster position by half a pixel for correct visualisation 
	// argh! wrong! correction has to happen in visualisation object, otherwise the actual data is wrong
	//header.x = header.x + (header.cellsize / 2);
	//header.y = header.y + (header.cellsize / 2);

	return true;
}

float* VtkRaster::loadDataFromASC(const std::string &fileName,
                                   double &x0,
                                   double &y0,
                                   size_t &width,
                                   size_t &height,
                                   double &delta)
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

		float* values = new float[header.ncols * header.nrows];

		int index(0);
		std::string s("");
		// read the file into a double-array
		for (int j = 0; j < header.nrows; j++)
		{
			index = (header.nrows - j - 1) * header.ncols;
			for (int i = 0; i < header.ncols; i++)
			{
				in >> s;
				values[index + i] = static_cast<float>(strtod(replaceString(",", ".", s).c_str(),0));
			}
		}

		in.close();
		return values;
	}
	return NULL;
}
/*
#ifdef libgeotiff_FOUND
bool VtkRaster::loadImageFromTIFF(const std::string &fileName,
                                  QImage &raster,
                                  QPointF &origin,
                                  double &cellsize)
{
	TIFF* tiff = XTIFFOpen(fileName.toStdString().c_str(), "r");

	if (tiff)
	{
		GTIF* geoTiff = GTIFNew(tiff);

		if (geoTiff)
		{
			int imgWidth = 0, imgHeight = 0, nImages = 0, pntCount = 0;
			double* pnts = 0;

			// get actual number of images in the tiff file
			do {
				nImages++;
			} while (TIFFReadDirectory(tiff));
			if (nImages > 1)
				std::cout << "OGSRaster::loadImageFromTIFF() - File contains " <<
				nImages << " images. This method is not tested for this case." <<
				std::endl;

			// get image size
			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH,  &imgWidth);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &imgHeight);

			// get cellsize
			// Note: GeoTiff allows anisotropic pixels. This is not supported here and equilateral pixels are assumed.
			if (TIFFGetField(tiff, GTIFF_PIXELSCALE, &pntCount, &pnts))
			{
				if (pnts[0] != pnts[1])
					std::cout <<
					"OGSRaster::loadImageFromTIFF() - Warning: Original raster data has anisotrop pixel size!"
					          << std::endl;
				cellsize = pnts[0];
			}

			// get upper left point / origin
			if (TIFFGetField(tiff, GTIFF_TIEPOINTS, &pntCount, &pnts))
			{
				origin.setX(pnts[3]);
				origin.setY(pnts[4] - (imgHeight * cellsize)); // the origin should be the lower left corner of the img
			}

			// read pixel values
			uint32* pixVal =
			        (uint32*) _TIFFmalloc(imgWidth * imgHeight * sizeof (uint32));
			if ((imgWidth > 0) && (imgHeight > 0))
				if (!TIFFReadRGBAImage(tiff, imgWidth, imgHeight, pixVal, 0))
				{
					std::cout <<
					"OGSRaster::loadImageFromTIFF() - Error reading GeoTIFF file."
					          << std::endl;
					_TIFFfree(pixVal);
					GTIFFree(geoTiff);
					XTIFFClose(tiff);
					return false;
				}

			// read colormap if it exists
			uint16* cmap_red = NULL, * cmap_green = NULL, * cmap_blue = NULL;
			int colormap_used = TIFFGetField(tiff,
			                                 TIFFTAG_COLORMAP,
			                                 &cmap_red,
			                                 &cmap_green,
			                                 &cmap_blue);

			int lineindex = 0, idx = 0;
			QImage img(imgWidth, imgHeight, QImage::Format_ARGB32);

			int* pxl (new int[4]);
			for (int j = 0; j < imgHeight; j++)
			{
				lineindex = j * imgWidth;
				for (int i = 0; i < imgWidth; i++)
				{ // scale intensities and set nodata values to white (i.e. the background colour)
					idx = TIFFGetR(pixVal[lineindex + i]);
					if (colormap_used)
						img.setPixel(i,j,
						             qRgba(cmap_red[idx] >> 8,
						                   cmap_green[idx] >> 8,
						                   cmap_blue[idx] >> 8, 255));
					else
					{
						//img.setPixel(i,j, qRgba(TIFFGetB(pixVal[idx]), TIFFGetG(pixVal[idx]), TIFFGetR(pixVal[idx]), TIFFGetA(pixVal[idx])));
						uint32toRGBA(pixVal[lineindex + i], pxl);
						img.setPixel(i,j,
						             qRgba(pxl[0], pxl[1], pxl[2], pxl[3]));
					}
				}
			}
			delete [] pxl;

			raster = img;

			_TIFFfree(pixVal);
			GTIFFree(geoTiff);
			XTIFFClose(tiff);
			return true;
		}

		XTIFFClose(tiff);
		std::cout <<
		"OGSRaster::loadImageFromTIFF() - File not recognised as GeoTIFF-Image." <<
		std::endl;
		return false;
	}

	std::cout << "OGSRaster::loadImageFromTIFF() - File not recognised as TIFF-Image." <<
	std::endl;
	return false;
}
#endif

bool VtkRaster::loadImageFromFile(const std::string &fileName, QImage &raster)
{
	return raster.load(fileName);
}

void VtkRaster::convertToGreyscale(vtkImageData &raster, const int &min, const int &max)
{
	int value = 0;
	double scalingFactor = 255.0 / (max - min);

	for (int i = 0; i < raster.width(); i++)
		for (int j = 0; j < raster.height(); j++)
		{
			QRgb pix = raster.pixel(i,j);
			value =
			        static_cast<int>(floor(((0.3 * qRed(pix) + 0.6 * qGreen(pix) +
			                                 0.1 *
			                                 qBlue(pix)) - min) * scalingFactor));
			raster.setPixel(i, j, qRgb(value, value, value));
		}
}

int* VtkRaster::getGreyscaleData(vtkImageData &raster, const int &min, const int &max)
{
	int index = 0;
	double scalingFactor = 255.0 / (max - min);
	int* pixVal (new int[raster.height() * raster.width()]);

	for (int j = 0; j < raster.height(); j++)
	{
		index = j * raster.width();
		for (int i = 0; i < raster.width(); i++)
		{
			QRgb pix = raster.pixel(i,j);
			pixVal[index +
			       i] =
			        static_cast<int>(floor(((0.3 * qRed(pix) + 0.6 * qGreen(pix) +
			                                 0.1 *
			                                 qBlue(pix)) - min) * scalingFactor));
		}
	}
	return pixVal;
}
*/

void VtkRaster::uint32toRGBA(const unsigned int s, int* p)
{
	p[3]   = s / (256 * 256 * 256);
	int r  = s % (256 * 256 * 256);
	p[2]   = r / (256 * 256);
	r     %= (256 * 256);
	p[1]   = r / 256;
	p[0]   = r % 256;
}
