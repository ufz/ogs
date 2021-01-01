/**
 * \file
 * \author Karsten Rink
 * \date   2012-02-01
 * \brief  Definition of the VtkRaster class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <string>
#include "GeoLib/Raster.h"

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
     * \brief Returns a VtkImageAlgorithm from an array of pixel values and some image meta data.
     */
    static vtkImageImport* loadImageFromArray(double const*const data_array, GeoLib::RasterHeader header);
    /**
     * \brief Loads an image- or raster-file into an vtkImageAlgorithm-Object.
     *
     * Public method for loading all data formats. Internally the method automatically differentiates between
     * images and georeferenced files and then calls the appropriate method for reading the file.
     * \param fileName Filename of the file that should be loaded.
     * \return The ImageAlgorithm-object.
     */
    static vtkImageAlgorithm* loadImage(const std::string &fileName);
private:
#ifdef GEOTIFF_FOUND
    /**
     * Loads ArcGIS asc-files to a QPixmap object and automatically does a contrast stretching to adjust values to 8 bit greyscale images.
     * \param fileName Filename of the file that should be loaded.
     * \return A vtkImageImport-object (derived from vtkImageAlgorithm).
     */
    static vtkImageAlgorithm* loadImageFromTIFF(const std::string& fileName);
#endif

    /**
     * Loads image files into a QPixmap object. Since images are not geo-referenced no origin point will be returned.
     * \param fileName Filename of the file that should be loaded.
     * \return vtkImageReader2-object containing the image data.
     */
    static vtkImageReader2* loadImageFromFile(const std::string &fileName);

    /**
     * Tries to find a world file for the image given by the filename.
     * World files can have a number of extensions depending on the programme
     * used to write the image and this method just cycles through the
     * possibilities, returning the first match it finds.
     */
    static std::string findWorldFile(const std::string& filename);

    /**
     * Tries to find and load the world file associated with a
     * BMP/JPG/PNG-file and create a RasterHeader for the image.
     */
    static bool readWorldFile(std::string const& filename, vtkImageReader2* image);
};
