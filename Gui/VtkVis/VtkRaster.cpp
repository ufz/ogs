/**
 * \file
 * \author Karsten Rink
 * \date   2012-02-01
 * \brief  Implementation of the VtkRaster class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkRaster.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <QFileInfo>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "StringTools.h"

// GeoLib
#include "Raster.h"

#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageReader2.h>
#include <vtkPNGReader.h>
#include <vtkJPEGReader.h>
#include <vtkBMPReader.h>

#ifdef libgeotiff_FOUND
#include "geo_tiffp.h"
#include "xtiffio.h"
#endif

vtkImageAlgorithm* VtkRaster::loadImage(const std::string &fileName,
                                        double& x0, double& y0, double& delta)
{
	QFileInfo fileInfo(QString::fromStdString(fileName));


	GeoLib::Raster *raster(nullptr);
	if (fileInfo.suffix().toLower() == "asc") {
		raster = GeoLib::Raster::getRasterFromASCFile(fileName);
	}
	else if (fileInfo.suffix().toLower() == "grd")
	{
		raster = GeoLib::Raster::getRasterFromSurferFile(fileName);
	}
	if (raster) {
		x0 = raster->getOrigin()[0];
		y0 = raster->getOrigin()[1];
		double const*const data (raster->begin());
		return VtkRaster::loadImageFromArray(data, x0, y0,
						raster->getNCols(), raster->getNRows(), raster->getRasterPixelDistance(),
						raster->getNoDataValue());
	}
#ifdef libgeotiff_FOUND
	else if ((fileInfo.suffix().toLower() == "tif") || (fileInfo.suffix().toLower() == "tiff"))
		return loadImageFromTIFF(fileName, x0, y0, delta);
#endif
	else
		return loadImageFromFile(fileName);
}

vtkImageImport* VtkRaster::loadImageFromArray(double const*const data_array, double x0, double y0, std::size_t width, std::size_t height, double delta, double noData)
{
	const unsigned length = height*width;
	float* data = new float[length*2];
	float max_val=noData;
	for (unsigned j=0; j<length; ++j)
	{
		data[j*2] = static_cast<float>(data_array[j]);
		max_val = (data[j*2]>max_val) ? data[j*2] : max_val;
	}
	for (unsigned j=0; j<length; ++j)
	{
		if (data[j*2]==noData)
		{
			data[j*2] = max_val;
			data[j*2+1] = 0;
		}
		else
		{
			//data[j*2] = max_val-data[j*2];//delete;
			data[j*2+1] = max_val;
		}
	}

	vtkImageImport* image = vtkImageImport::New();
		image->SetDataSpacing(delta, delta, delta);
		image->SetDataOrigin(x0+(delta/2.0), y0+(delta/2.0), 0);	// translate whole mesh by half a pixel in x and y
		image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
		image->SetDataExtent(0, width-1, 0, height-1-1, 0, 0);
		image->SetDataExtentToWholeExtent();
		image->SetDataScalarTypeToFloat();
		image->SetNumberOfScalarComponents(2);
		image->SetImportVoidPointer(data, 0);
		image->Update();

	return image;
}




#ifdef libgeotiff_FOUND
vtkImageImport* VtkRaster::loadImageFromTIFF(const std::string &fileName,
                                  double &x0, double &y0,
                                  double &cellsize)
{
	TIFF* tiff = XTIFFOpen(fileName.c_str(), "r");

	if (tiff)
	{
		GTIF* geoTiff = GTIFNew(tiff);

		if (geoTiff)
		{
			int imgWidth = 0, imgHeight = 0, nImages = 0, pntCount = 0;
			double* pnts = 0;

			// get actual number of images in the tiff file
			do {
				++nImages;
			} while (TIFFReadDirectory(tiff));
			if (nImages > 1)
				INFO("VtkRaster::loadImageFromTIFF() - File contains %d images. This method is not tested for this case.", nImages);

			// get image size
			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH,  &imgWidth);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &imgHeight);

			// get cellsize
			// Note: GeoTiff allows anisotropic pixels. This is not supported here and equilateral pixels are assumed.
			if (TIFFGetField(tiff, GTIFF_PIXELSCALE, &pntCount, &pnts))
			{
				if (pnts[0] != pnts[1])
					WARN("VtkRaster::loadImageFromTIFF(): Original raster data has anisotrop pixel size!");
				cellsize = pnts[0];
			}

			// get upper left point / origin
			if (TIFFGetField(tiff, GTIFF_TIEPOINTS, &pntCount, &pnts))
			{
				x0 = pnts[3];
				y0 = pnts[4] - (imgHeight * cellsize); // the origin should be the lower left corner of the img
			}

			// read pixel values
			uint32* pixVal =
			        (uint32*) _TIFFmalloc(imgWidth * imgHeight * sizeof (uint32));
			if ((imgWidth > 0) && (imgHeight > 0))
				if (!TIFFReadRGBAImage(tiff, imgWidth, imgHeight, pixVal, 0))
				{
					ERR("VtkRaster::loadImageFromTIFF(): reading GeoTIFF file.");
					_TIFFfree(pixVal);
					GTIFFree(geoTiff);
					XTIFFClose(tiff);
					return nullptr;
				}

			// check for colormap
			uint16 photometric;
			TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);
			// read colormap
			uint16* cmap_red = nullptr, * cmap_green = nullptr, * cmap_blue = nullptr;
			int colormap_used = TIFFGetField(tiff,
											TIFFTAG_COLORMAP,
											&cmap_red,
											&cmap_green,
											&cmap_blue);

			float* data = new float[imgWidth * imgHeight * 4];
			int* pxl (new int[4]);
			for (int j = 0; j < imgHeight; ++j)
			{
				int lineindex = j * imgWidth;
				for (int i = 0; i < imgWidth; ++i)
				{ // scale intensities and set nodata values to white (i.e. the background colour)
					unsigned pxl_idx(lineindex+i);
					unsigned pos  = 4 * (pxl_idx);
					if (photometric==1 && colormap_used==1)
					{
						int idx = TIFFGetR(pixVal[pxl_idx]);
						data[pos]   = cmap_red[idx] >> 8;
						data[pos+1] = cmap_green[idx] >> 8;
						data[pos+2] = cmap_blue[idx] >> 8;
						data[pos+3] = 1;
					}
					else
					{
						data[pos]   = static_cast<float>(TIFFGetR(pixVal[pxl_idx]));
						data[pos+1] = static_cast<float>(TIFFGetG(pixVal[pxl_idx]));
						data[pos+2] = static_cast<float>(TIFFGetB(pixVal[pxl_idx]));
						data[pos+3] = static_cast<float>(TIFFGetA(pixVal[pxl_idx]));
					}
				}
			}
			delete [] pxl;

			// set transparency values according to maximum pixel value
			if (photometric==1)
			{
				float max_val(0);
				unsigned nPixels = 4*imgWidth*imgHeight;
				for (unsigned j = 0; j < nPixels; ++j)
					if (data[j]>max_val)
						max_val = data[j];

				for (unsigned j = 0; j < nPixels; j+=4)
					data[j+3] = max_val;
			}

			vtkImageImport* image = vtkImageImport::New();
				image->SetDataOrigin(x0, y0, 0);
				image->SetDataSpacing(cellsize, cellsize, cellsize);
				image->SetWholeExtent(0, imgWidth-1, 0, imgHeight-1, 0, 0);
				image->SetDataExtent(0, imgWidth-1, 0, imgHeight-1, 0, 0);
				image->SetDataExtentToWholeExtent();
				image->SetDataScalarTypeToFloat();
				image->SetNumberOfScalarComponents(4);
				image->SetImportVoidPointer(data, 0);
				image->Update();

			_TIFFfree(pixVal);
			GTIFFree(geoTiff);
			XTIFFClose(tiff);
			return image;
		}

		XTIFFClose(tiff);
		ERR("VtkRaster::loadImageFromTIFF() - File not recognised as GeoTIFF-Image.")
		return nullptr;
	}

	ERR("VtkRaster::loadImageFromTIFF() - File not recognised as TIFF-Image.");
	return nullptr;
}
#endif

vtkImageReader2* VtkRaster::loadImageFromFile(const std::string &fileName)
{
	QString file_name (QString::fromStdString(fileName));
	QFileInfo fi(file_name);
	vtkImageReader2* image(nullptr);

	if (fi.suffix().toLower() == "png")
		image = vtkPNGReader::New();
	else if ((fi.suffix().toLower() == "jpg") || (fi.suffix().toLower() == "jpeg"))
		image = vtkJPEGReader::New();
	else if (fi.suffix().toLower() == "bmp")
		image = vtkBMPReader::New();
	else
	{
		ERR("VtkRaster::readImageFromFile(): File format not support, please convert to BMP, JPG, PNG or TIFF.");
		return nullptr;
	}

	image->SetFileName(fileName.c_str());
	image->GetOutput()->AllocateScalars(VTK_FLOAT, 1);
	image->Update();
	return image;
}

void VtkRaster::uint32toRGBA(const unsigned int s, int* p)
{
	p[3]   = s / (256 * 256 * 256);
	int r  = s % (256 * 256 * 256);
	p[2]   = r / (256 * 256);
	r     %= (256 * 256);
	p[1]   = r / 256;
	p[0]   = r % 256;
}
