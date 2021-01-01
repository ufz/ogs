/**
 * \file
 * \author Karsten Rink
 * \date   2012-02-01
 * \brief  Implementation of the VtkRaster class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkRaster.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <QFileInfo>

#include <vtkBMPReader.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageReader2.h>
#include <vtkJPEGReader.h>
#include <vtkPNGReader.h>
#include <vtkTIFFReader.h>

#ifdef GEOTIFF_FOUND
#include "geo_tiffp.h"
#include "xtiffio.h"
#endif

#include <memory>

#include "Applications/FileIO/AsciiRasterInterface.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/Raster.h"

vtkImageAlgorithm* VtkRaster::loadImage(const std::string &fileName)
{
    QFileInfo fileInfo(QString::fromStdString(fileName));

    std::unique_ptr<GeoLib::Raster> raster(nullptr);
    if (fileInfo.suffix().toLower() == "asc")
    {
        raster.reset(FileIO::AsciiRasterInterface::getRasterFromASCFile(fileName));
    }
    else if (fileInfo.suffix().toLower() == "grd")
    {
        raster.reset(
            FileIO::AsciiRasterInterface::getRasterFromSurferFile(fileName));
    }
    if (raster)
    {
        return VtkRaster::loadImageFromArray(raster->begin(),
                                             raster->getHeader());
    }
    if ((fileInfo.suffix().toLower() == "tif") ||
        (fileInfo.suffix().toLower() == "tiff"))
    {
#ifdef GEOTIFF_FOUND
        return loadImageFromTIFF(fileName);
#else
        ERR("VtkRaster::loadImage(): GeoTiff file format not supported in this "
            "version! Trying to parse as Tiff-file.");
        return loadImageFromFile(fileName);
#endif
    }

    return loadImageFromFile(fileName);
}
vtkImageImport* VtkRaster::loadImageFromArray(double const*const data_array, GeoLib::RasterHeader header)
{
    const unsigned length = header.n_rows * header.n_cols * header.n_depth;
    auto* data = new float[length * 2];
    float max_val = static_cast<float>(*std::max_element(data_array, data_array+length));
    for (unsigned j=0; j<length; ++j)
    {
        data[j*2] = static_cast<float>(data_array[j]);
        if (fabs(data[j*2]-header.no_data) < std::numeric_limits<double>::epsilon())
        {
            data[j*2] = max_val;
            data[j*2+1] = 0;
        }
        else
        {
            data[j * 2 + 1] = max_val;
        }
    }

    vtkImageImport* image = vtkImageImport::New();
        image->SetDataSpacing(header.cell_size, header.cell_size, header.cell_size);
        image->SetDataOrigin(header.origin[0]+(header.cell_size/2.0), header.origin[1]+(header.cell_size/2.0), 0);    // translate whole mesh by half a pixel in x and y
        image->SetWholeExtent(0, header.n_cols - 1, 0, header.n_rows - 1, 0,
                              header.n_depth - 1);
        image->SetDataExtent(0, header.n_cols - 1, 0, header.n_rows - 1, 0,
                             header.n_depth - 1);
        image->SetDataExtentToWholeExtent();
        image->SetDataScalarTypeToFloat();
        image->SetNumberOfScalarComponents(2);
        image->SetImportVoidPointer(data, 0);
        image->Update();

    return image;
}

#ifdef GEOTIFF_FOUND
vtkImageAlgorithm* VtkRaster::loadImageFromTIFF(const std::string& fileName)
{
    TIFF* tiff = XTIFFOpen(fileName.c_str(), "r");

    if (tiff)
    {
        GTIF* geoTiff = GTIFNew(tiff);

        int version[3];
        int count (0);
        GTIFDirectoryInfo(geoTiff, version, &count);

        if (count == 0)
            WARN("VtkRaster::loadImageFromTIFF - file is not georeferenced.");

        if (geoTiff)
        {
            double x0 = 0.0;
            double y0 = 0.0;
            double cellsize = 1.0;
            int imgWidth = 0;
            int imgHeight = 0;
            int nImages = 0;
            int pntCount = 0;
            double* pnts = nullptr;

            // get actual number of images in the tiff file
            do {
                ++nImages;
            } while (TIFFReadDirectory(tiff));
            if (nImages > 1)
                INFO(
                    "VtkRaster::loadImageFromTIFF() - File contains {:d} "
                    "images. This method is not tested for this case.",
                    nImages);

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
            auto* pixVal = static_cast<uint32*>(
                _TIFFmalloc(imgWidth * imgHeight * sizeof(uint32)));
            if ((imgWidth > 0) && (imgHeight > 0))
            {
                if (!TIFFReadRGBAImage(tiff, imgWidth, imgHeight, pixVal, 0))
                {
                    ERR("VtkRaster::loadImageFromTIFF(): reading GeoTIFF file.");
                    _TIFFfree(pixVal);
                    GTIFFree(geoTiff);
                    XTIFFClose(tiff);
                    return nullptr;
                }
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

            auto* data = new float[imgWidth * imgHeight * 4];
            auto* pxl(new int[4]);
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
                        data[pos]   = static_cast<float>(cmap_red[idx] >> 8);
                        data[pos+1] = static_cast<float>(cmap_green[idx] >> 8);
                        data[pos+2] = static_cast<float>(cmap_blue[idx] >> 8);
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
                {
                    if (data[j] > max_val)
                    {
                        max_val = data[j];
                    }
                }

                for (unsigned j = 0; j < nPixels; j += 4)
                {
                    data[j + 3] = max_val;
                }
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
        ERR("VtkRaster::loadImageFromTIFF() - File not recognised as "
            "GeoTIFF-Image.");
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
    {
        image = vtkPNGReader::New();
    }
    else if ((fi.suffix().toLower() == "tif") ||
             (fi.suffix().toLower() == "tiff"))
    {
        image = vtkTIFFReader::New();
    }
    else if ((fi.suffix().toLower() == "jpg") ||
             (fi.suffix().toLower() == "jpeg"))
    {
        image = vtkJPEGReader::New();
    }
    else if (fi.suffix().toLower() == "bmp")
    {
        image = vtkBMPReader::New();
    }
    else
    {
        ERR("VtkRaster::readImageFromFile(): File format not supported, please "
            "convert to BMP, JPG, PNG or TIFF.");
        return nullptr;
    }

    image->SetFileName(fileName.c_str());
    image->GetOutput()->AllocateScalars(VTK_FLOAT, 1);
    image->Update();
    readWorldFile(fileName, image);
    return image;
}

std::string VtkRaster::findWorldFile(std::string const& filename)
{
    std::string const no_ext = BaseLib::dropFileExtension(filename);

    constexpr std::array supported_extensions = {
        ".pgw",  ".pngw", ".pgwx", ".jgw",  ".jpgw", ".jgwx", ".tfw",
        ".tifw", ".tfwx", ".bpw",  ".bmpw", ".bpwx", ".wld"};

    auto const res =
        std::find_if(supported_extensions.begin(), supported_extensions.end(),
                     [&no_ext](auto const& ext) -> bool {
                         return BaseLib::IsFileExisting(no_ext + ext);
                     });
    if (res != supported_extensions.end())
    {
        return no_ext + *res;
    }

    // no world file found
    return {};
}

bool VtkRaster::readWorldFile(std::string const& filename,
                              vtkImageReader2* image)
{
    std::string const world_file = findWorldFile(filename);
    if (world_file.empty())
    {
        WARN("No world file found. Image is not georeferenced.");
        return false;
    }

    std::ifstream in(world_file.c_str());
    if (!in.is_open())
    {
        ERR("VtkRaster::readWorldFile(): Could not open file {:s}.", filename);
        return false;
    }

    std::string line;
    // x-scaling
    if (!std::getline(in, line))
    {
        return false;
    }
    double const delta_x = BaseLib::str2number<double>(line);
    // 2x rotation (ignored)
    if (!(std::getline(in, line) && std::getline(in, line)))
    {
        return false;
    }
    // negative y-scaling
    if (!std::getline(in, line))
    {
        return false;
    }
    double const delta_y = BaseLib::str2number<double>(line);
    if (delta_x != -delta_y)
        WARN(
            "Anisotropic pixel size detected ({:f} vs {:f}). An isotropic "
            "spacing of {:f} is assumed, be aware results may be wrong.",
            delta_x, delta_y, delta_x);
    // x-translation
    if (!std::getline(in, line))
    {
        return false;
    }
    double const x0 = BaseLib::str2number<double>(line);
    // y-translation
    if (!std::getline(in, line))
    {
        return false;
    }
    double const y0 = BaseLib::str2number<double>(line);

    int extent[3];
    image->GetOutput()->GetDimensions(extent);
    image->SetDataSpacing(delta_x, delta_x, delta_x);
    // for GIS the origin is on the lower left, for VTK it is on the upper left
    double const vtk_y0 = y0 + (extent[1] * delta_y);
    image->SetDataOrigin(x0, vtk_y0, 0);
    return true;
}
