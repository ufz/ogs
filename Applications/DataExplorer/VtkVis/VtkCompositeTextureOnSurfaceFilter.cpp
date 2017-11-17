/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-18
 * \brief  Implementation of the VtkCompositeTextureOnSurfaceFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeTextureOnSurfaceFilter.h"

#include <logog/include/logog.hpp>

#include "VtkTextureOnSurfaceFilter.h"
#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include "VtkGeoImageSource.h"
#include "VtkRaster.h"
#include "NetCdfDialog/NetCdfConfigureDialog.h"

#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>

//#include "VtkCompositeColormapToImageFilter.h"


VtkCompositeTextureOnSurfaceFilter::VtkCompositeTextureOnSurfaceFilter(
        vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm)
{
    this->init();
}

VtkCompositeTextureOnSurfaceFilter::~VtkCompositeTextureOnSurfaceFilter() =
    default;

void VtkCompositeTextureOnSurfaceFilter::init()
{
    this->_inputDataObjectType = VTK_DATA_SET;
    this->_outputDataObjectType = VTK_POLY_DATA;

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter;
    VtkTextureOnSurfaceFilter* surface = VtkTextureOnSurfaceFilter::New();

    if (dynamic_cast<vtkUnstructuredGrid*>(_inputAlgorithm->GetOutputDataObject(0)))
    {
        surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfaceFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
        surface->SetInputConnection(surfaceFilter->GetOutputPort());
    }
    else
        surface->SetInputConnection(_inputAlgorithm->GetOutputPort());

    QWidget* parent = nullptr;
    QSettings settings;
    QString fileName = QFileDialog::getOpenFileName(parent, "Select raster file to apply as texture",
                                                    settings.value("lastOpenedTextureFileDirectory").toString(),
                                                    "Raster files (*.asc *.grd *.bmp *.jpg *.png *.tif);;NetCDF files (*.nc);;");
    QFileInfo fi(fileName);

    if ((fi.suffix().toLower() == "asc") || (fi.suffix().toLower() == "tif") ||
        (fi.suffix().toLower() == "png") || (fi.suffix().toLower() == "grd") ||
        (fi.suffix().toLower() == "jpg") || (fi.suffix().toLower() == "bmp"))
    {
        double x0(0), y0(0), scalingFactor(1);
        std::string name = fileName.toStdString();
        vtkImageAlgorithm* image = VtkRaster::loadImage(name, x0, y0, scalingFactor);
        surface->SetRaster(image, x0, y0, scalingFactor);
        surface->Update();

        QDir dir = QDir(fileName);
        settings.setValue("lastOpenedTextureFileDirectory", dir.absolutePath());
    }
    else if (fi.suffix().toLower() == "nc")
    {
        NetCdfConfigureDialog dlg(fileName.toStdString().c_str());
        dlg.exec();
        if (dlg.getRaster() != nullptr)
        {
            VtkGeoImageSource* image = dlg.getRaster();
            double origin[3];
            image->GetOutput()->GetOrigin(origin);
            double spacing[3];
            image->GetOutput()->GetSpacing(spacing);
            surface->SetRaster(image, origin[0], origin[1], spacing[0]);
            surface->Update();
        }
    }
    else
        ERR("VtkCompositeTextureOnSurfaceFilter::init(): Error reading texture file.");

    _outputAlgorithm = surface;
}

void VtkCompositeTextureOnSurfaceFilter::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);
}
