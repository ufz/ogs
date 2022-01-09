/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-18
 * \brief  Implementation of the VtkCompositeTextureOnSurfaceFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeTextureOnSurfaceFilter.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "BaseLib/Logging.h"
#include "VtkGeoImageSource.h"
#include "VtkRaster.h"
#include "VtkTextureOnSurfaceFilter.h"
#ifdef OGS_USE_NETCDF
#include "NetCdfConfigureDialog.h"
#endif  // OGS_USE_NETCDF

#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>

//#include "VtkCompositeColormapToImageFilter.h"

VtkCompositeTextureOnSurfaceFilter::VtkCompositeTextureOnSurfaceFilter(
    vtkAlgorithm* inputAlgorithm)
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

    if (dynamic_cast<vtkUnstructuredGrid*>(
            _inputAlgorithm->GetOutputDataObject(0)))
    {
        surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfaceFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
        surface->SetInputConnection(surfaceFilter->GetOutputPort());
    }
    else
    {
        surface->SetInputConnection(_inputAlgorithm->GetOutputPort());
    }

    QWidget* parent = nullptr;
    QSettings settings;
    QString fileName = QFileDialog::getOpenFileName(
        parent, "Select raster file to apply as texture",
        settings.value("lastOpenedTextureFileDirectory").toString(),
        "Raster files (*.asc *.grd *.bmp *.jpg *.png *.tif);;"
#ifdef OGS_USE_NETCDF
        "NetCDF files (*.nc);;"
#endif  // OGS_USE_NETCDF
    );
    QFileInfo fi(fileName);

    if ((fi.suffix().toLower() == "asc") || (fi.suffix().toLower() == "tif") ||
        (fi.suffix().toLower() == "png") || (fi.suffix().toLower() == "grd") ||
        (fi.suffix().toLower() == "jpg") || (fi.suffix().toLower() == "bmp"))
    {
        std::string name = fileName.toStdString();
        vtkImageAlgorithm* image = VtkRaster::loadImage(name);
        surface->SetRaster(image);
        surface->Update();

        QDir dir = QDir(fileName);
        settings.setValue("lastOpenedTextureFileDirectory", dir.absolutePath());
    }
#ifdef OGS_USE_NETCDF
    else if (fi.suffix().toLower() == "nc")
    {
        NetCdfConfigureDialog dlg(fileName.toStdString().c_str());
        dlg.exec();
        if (dlg.getRaster() != nullptr)
        {
            VtkGeoImageSource* image = dlg.getRaster();
            surface->SetRaster(image);
            surface->Update();
        }
    }
#endif  // OGS_USE_NETCDF
    else
        ERR("VtkCompositeTextureOnSurfaceFilter::init(): Error reading texture "
            "file.");

    _outputAlgorithm = surface;
}

void VtkCompositeTextureOnSurfaceFilter::SetUserProperty(QString name,
                                                         QVariant value)
{
    VtkAlgorithmProperties::SetUserProperty(name, value);
}
