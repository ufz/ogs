/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VtkRaster.cpp
 *
 * Created on 2012-02-01 by Karsten Rink
 */

#include "VtkRaster.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <QFileInfo>

#include "StringTools.h"

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

	if (fileInfo.suffix().toLower() == "asc" || fileInfo.suffix().toLower() == "grd")
        return loadImageFromASC(fileName, x0, y0, delta);
#ifdef libgeotiff_FOUND
	else if ((fileInfo.suffix().toLower() == "tif") || (fileInfo.suffix().toLower() == "tiff"))
		return loadImageFromTIFF(fileName, x0, y0, delta);
#endif
	else
		return loadImageFromFile(fileName);
}

vtkImageImport* VtkRaster::loadImageFromASC(const std::string &fileName,
                                            double& x0, double& y0, double& delta)
{
	unsigned width(0), height(0);
	float* data;

	if (fileName.substr(fileName.length()-3, 3).compare("asc") == 0)
		data = loadDataFromASC(fileName, x0, y0, width, height, delta);
	else
		data = loadDataFromSurfer(fileName, x0, y0, width, height, delta);

	vtkImageImport* image = vtkImageImport::New();
		image->SetDataSpacing(delta, delta,delta);
		image->SetDataOrigin(x0+(delta/2.0), y0+(delta/2.0), 0);	// translate whole mesh by half a pixel in x and y
		image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
		image->SetDataExtent(0, width-1, 0, height-1, 0, 0);
		image->SetDataExtentToWholeExtent();
		image->SetDataScalarTypeToFloat();
		image->SetNumberOfScalarComponents(2);
		image->SetImportVoidPointer(data, 0);
		image->Update();

	return image;
}

vtkImageImport* VtkRaster::loadImageFromArray(double* data_array, double &x0, double &y0, unsigned &width, unsigned &height, double &delta, double noData)
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
		header.noData = value.c_str();
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
                                   unsigned &width,
                                   unsigned &height,
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

		float* values = new float[header.ncols * header.nrows * 2];

		int col_index(0);
		int noData = atoi(header.noData.c_str());
		float max_val = noData;
		std::string s("");
		// read the file into a double-array
		for (int j = 0; j < header.nrows; ++j)
		{
			col_index = (header.nrows - j - 1) * header.ncols;
			for (int i = 0; i < header.ncols; ++i)
			{
				in >> s;
				unsigned index = 2*(col_index+i);
				values[index] = static_cast<float>(strtod(replaceString(",", ".", s).c_str(),0));
				if (values[index] > max_val)
					max_val = values[index];
			}
		}

		// shift noData values into normal pixel-range and set transparancy values for all pixels
		unsigned nPixels = header.ncols * header.nrows;
		for (unsigned j = 0; j < nPixels; ++j)
		{
			if (values[j*2] == noData)
			{
				values[j*2] = max_val;
				values[j*2+1] = 0;
			}
			else
				values[j*2+1] = max_val;
		}

		in.close();
		return values;
	}
	return NULL;
}

bool VtkRaster::readSurferHeader(ascHeader &header, std::ifstream &in)
{
	std::string line, tag, value;
	double min, max;

	in >> tag;

	if (tag.compare("DSAA") != 0)
	{
		std::cout << "Error in readSurferHeader() - No Surfer file..." << std::endl;
		return false;
	}
	else
	{
		in >> header.ncols >> header.nrows;
		in >> min >> max;
		header.x = min;
		header.cellsize = (max-min)/(double)header.ncols;

		in >> min >> max;
		header.y = min;

		if (ceil((max-min)/(double)header.nrows) == ceil(header.cellsize))
			header.cellsize = ceil(header.cellsize);
		else
		{
			std::cout << "Error in readSurferHeader() - Anisotropic cellsize detected..." << std::endl;
			return 0;
		}
		in >> min >> max; // ignore min- and max-values

		header.noData = "1.70141E+038";
	}

	return true;
}

float* VtkRaster::loadDataFromSurfer(const std::string &fileName,
                                   double &x0,
                                   double &y0,
                                   unsigned &width,
                                   unsigned &height,
                                   double &delta)
{
	std::ifstream in( fileName.c_str() );

	if (!in.is_open())
	{
		std::cout << "VtkRaster::loadImageFromSurfer() - Could not open file..." << std::endl;
		return NULL;
	}

	ascHeader header;

	if (readSurferHeader(header, in))
	{
		x0     = header.x;
		y0     = header.y;
		width  = header.ncols;
		height = header.nrows;
		delta  = header.cellsize;

		float* values = new float[header.ncols * header.nrows * 2];

		int col_index(0);
		int noData = -9999;
		float max_val = noData;
		std::string s("");
		// read the file into a double-array
		for (int j = 0; j < header.nrows; ++j)
		{
			col_index = j * header.ncols;
			for (int i = 0; i < header.ncols; ++i)
			{
				in >> s;
				if (s.compare(header.noData) == 0)
					s = "-9999";
				unsigned index = 2*(col_index+i);
				values[index] = static_cast<float>(strtod(replaceString(",", ".", s).c_str(),0));
				if (values[index] > max_val)
					max_val = values[index];
			}
		}

		// shift noData values into normal pixel-range and set transparancy values for all pixels
		unsigned nPixels = header.ncols * header.nrows;
		for (unsigned j = 0; j < nPixels; ++j)
		{
			if (values[j*2] == noData)
			{
				values[j*2] = max_val;
				values[j*2+1] = 0;
			}
			else
				values[j*2+1] = max_val;
		}

		in.close();
		return values;
	}
	return NULL;
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
				std::cout << "VtkRaster::loadImageFromTIFF() - File contains " <<
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
					"VtkRaster::loadImageFromTIFF() - Warning: Original raster data has anisotrop pixel size!"
					          << std::endl;
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
					std::cout <<
					"VtkRaster::loadImageFromTIFF() - Error reading GeoTIFF file."
					          << std::endl;
					_TIFFfree(pixVal);
					GTIFFree(geoTiff);
					XTIFFClose(tiff);
					return NULL;
				}

			// check for colormap
			uint16 photometric;
			TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);
			// read colormap
			uint16* cmap_red = NULL, * cmap_green = NULL, * cmap_blue = NULL;
			/*int colormap_used = */TIFFGetField(tiff,
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
					if (photometric==1)
					{
						int idx = TIFFGetR(pixVal[pxl_idx]);
						data[pos]   = cmap_red[idx] >> 8;
						data[pos+1] = cmap_green[idx] >> 8;
						data[pos+2] = cmap_blue[idx] >> 8;
						data[pos+3] = 1;
					}
					else
					{
						data[pos]   = TIFFGetR(pixVal[pxl_idx]);
						data[pos+1] = TIFFGetG(pixVal[pxl_idx]);
						data[pos+2] = TIFFGetB(pixVal[pxl_idx]);
						data[pos+3] = TIFFGetA(pixVal[pxl_idx]);
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
		std::cout <<
		"VtkRaster::loadImageFromTIFF() - File not recognised as GeoTIFF-Image." <<
		std::endl;
		return NULL;
	}

	std::cout << "VtkRaster::loadImageFromTIFF() - File not recognised as TIFF-Image." <<
	std::endl;
	return NULL;
}
#endif

vtkImageReader2* VtkRaster::loadImageFromFile(const std::string &fileName)
{
	QString file_name (QString::fromStdString(fileName));
	QFileInfo fi(file_name);
	vtkImageReader2* image(NULL);

	if (fi.suffix().toLower() == "png")
		image = vtkPNGReader::New();
	else if ((fi.suffix().toLower() == "jpg") || (fi.suffix().toLower() == "jpeg"))
		image = vtkJPEGReader::New();
	else if (fi.suffix().toLower() == "bmp")
		image = vtkBMPReader::New();
	else
	{
		std::cout << "VtkRaster::readImageFromFile() - File format not support, please convert to BMP, JPG, PNG or TIFF..." << std::endl;
		return NULL;
	}

	image->SetFileName(fileName.c_str());
	image->GetOutput()->SetScalarTypeToFloat();
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
