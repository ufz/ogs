/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-20
 * \brief  Implementation of the VtkFilterFactory class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkFilterFactory.h"

#include "VtkCompositeColorByHeightFilter.h"
#include "VtkCompositeColormapToImageFilter.h"
#include "VtkCompositeContourFilter.h"
#include "VtkCompositeElementSelectionFilter.h"
#include "VtkCompositeGeoObjectFilter.h"
#include "VtkCompositeImageToCylindersFilter.h"
#include "VtkCompositeLineToTubeFilter.h"
#include "VtkCompositeNodeSelectionFilter.h"
#include "VtkCompositePointToGlyphFilter.h"
#include "VtkCompositeTextureOnSurfaceFilter.h"
#include "VtkCompositeThresholdFilter.h"
#include "VtkImageDataToLinePolyDataFilter.h"
#include "VtkCompositeImageToPointCloudFilter.h"

#include <vtkDataSetSurfaceFilter.h>

const QVector<VtkFilterInfo> VtkFilterFactory::GetFilterList()
{
    QVector<VtkFilterInfo> filterList;

    // Composite filters
    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeImageToCylindersFilter",
                                 "Image to bar chart",
                                 "This filter converts the red pixel values of the image into a bar graph.",
                                 VTK_IMAGE_DATA, VTK_POLY_DATA));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositePointToGlyphFilter",
                                 "Points to spheres",
                                 "This filter generates spheres on point data that can be scaled and colored by scalar data.",
                                 VTK_POLY_DATA, VTK_POLY_DATA));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeLineToTubeFilter",
                                 "Lines to tubes",
                                 "This filter will convert lines to tubes that can be colored by scalar data.",
                                 VTK_POLY_DATA, VTK_POLY_DATA));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeColormapToImageFilter",
                                 "Apply lookup table to image",
                                 "This filter will take an input image of any valid scalar type, and map the first component of the image through a lookup table.",
                                 VTK_IMAGE_DATA, VTK_IMAGE_DATA));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeTextureOnSurfaceFilter",
                                 "Apply texture to surface",
                                 "This filter assigns an image or raster file as a texture for the given surface.",
                                 VTK_POINT_SET, VTK_POLY_DATA));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeThresholdFilter",
                                 "Extract cells by threshold",
                                 "This filter extracts cells from any dataset type that satisfy a threshold criterion. A cell satisfies the criterion if the (first) scalar value of (every or any) point satisfies the criterion. For example this can be used to show only certain material groups in a mesh.",
                                 VTK_POINT_SET, VTK_UNSTRUCTURED_GRID));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeColorByHeightFilter",
                                 "Elevation-based colouring",
                                 "This filter will generate scalar values based on the elevation of each point in the dataset.",
                                 VTK_POINT_SET, VTK_POLY_DATA));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeContourFilter",
                                 "Generate contours based on scalar fields",
                                 "Visualisation of contour-lines/-planes within dense scalar fields.",
                                 VTK_UNSTRUCTURED_GRID, VTK_UNSTRUCTURED_GRID));

    filterList.push_back(VtkFilterInfo(
                                 "VtkCompositeImageToPointCloudFilter",
                                 "Image to point cloud",
                                 "This filter creates point clouds with a density based on the first component of pixel values.",
                                 VTK_IMAGE_DATA, VTK_POLY_DATA));

    // Simple filters
    filterList.push_back(VtkFilterInfo(
                                 "VtkImageDataToLinePolyDataFilter",
                                 "Image to vertical lines",
                                 "This filter converts the red pixel values of the image to lines with length of the value.",
                                 VTK_IMAGE_DATA, VTK_POLY_DATA));

    // Standard VTK filter without properties
    filterList.push_back(VtkFilterInfo(
                                 "vtkDataSetSurfaceFilter",
                                 "Surface filter",
                                 "Extracts outer (polygonal) surface.",
                                 VTK_UNSTRUCTURED_GRID, VTK_POLY_DATA));

//      filterList.push_back(VtkFilterInfo(
//              "VtkCompositeSelectionFilter",
//              "Mesh Quality Filter",
//              "This filter calculates the quality of meshes and highlights deformed elements.",
//              VTK_UNSTRUCTURED_GRID, VTK_UNSTRUCTURED_GRID));

    return filterList;
}

VtkCompositeFilter* VtkFilterFactory::CreateCompositeFilter( QString type,
                                                             vtkAlgorithm* inputAlgorithm )
{
    if (type.compare(QString("VtkCompositeImageToCylindersFilter")) == 0)
        return new VtkCompositeImageToCylindersFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositePointToGlyphFilter")) == 0)
        return new VtkCompositePointToGlyphFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeLineToTubeFilter")) == 0)
        return new VtkCompositeLineToTubeFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeColormapToImageFilter")) == 0)
        return new VtkCompositeColormapToImageFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeTextureOnSurfaceFilter")) == 0)
        return new VtkCompositeTextureOnSurfaceFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeThresholdFilter")) == 0)
        return new VtkCompositeThresholdFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeColorByHeightFilter")) == 0)
        return new VtkCompositeColorByHeightFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeElementSelectionFilter")) == 0)
        return new VtkCompositeElementSelectionFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeNodeSelectionFilter")) == 0)
        return new VtkCompositeNodeSelectionFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeContourFilter")) == 0)
        return new VtkCompositeContourFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeGeoObjectFilter")) == 0)
        return new VtkCompositeGeoObjectFilter(inputAlgorithm);
    if (type.compare(QString("VtkCompositeImageToPointCloudFilter")) == 0)
        return new VtkCompositeImageToPointCloudFilter(inputAlgorithm);

    return nullptr;
}

vtkAlgorithm* VtkFilterFactory::CreateSimpleFilter( QString type )
{
    if (type.compare(QString("VtkImageDataToLinePolyDataFilter")) == 0)
        return VtkImageDataToLinePolyDataFilter::New();
    if (type.compare(QString("vtkDataSetSurfaceFilter")) == 0)
        return vtkDataSetSurfaceFilter::New();

    return nullptr;
}
