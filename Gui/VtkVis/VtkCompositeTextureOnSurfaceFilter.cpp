/**
 * \file VtkCompositeTextureOnSurfaceFilter.cpp
 * 18/11/2010 KR Initial implementation
 *
 * Implementation of VtkCompositeTextureOnSurfaceFilter class
 */

// ** INCLUDES **
#include "VtkCompositeTextureOnSurfaceFilter.h"
#include "VtkTextureOnSurfaceFilter.h"
#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include "VtkGeoImageSource.h"
#include "VtkRaster.h"
#include "NetCdfConfigureDialog.h"

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

VtkCompositeTextureOnSurfaceFilter::~VtkCompositeTextureOnSurfaceFilter()
{
}

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

	QWidget* parent = 0;
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getOpenFileName(parent,
	                                                "Select raster file to apply as texture",
	                                                settings.value("lastOpenedTextureFileDirectory").
	                                                toString(),
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
		if (dlg.getRaster() != NULL)
		{
			VtkGeoImageSource* image = dlg.getRaster();
			double origin[3];
			image->GetOutput()->GetOrigin(origin);
			double spacing[3];
			image->GetOutput()->GetSpacing(spacing);
/*
			VtkCompositeColormapToImageFilter* cm = new VtkCompositeColormapToImageFilter(image);
			vtkImageAlgorithm* img = dynamic_cast<vtkImageAlgorithm*>(cm->GetOutputAlgorithm());
*/
			surface->SetRaster(image, origin[0], origin[1], spacing[0]);
			surface->Update();
		}
	}
	else
		std::cout <<
		"VtkCompositeTextureOnSurfaceFilter.init() - Error reading texture file..." <<
		std::endl;

	_outputAlgorithm = surface;
}

void VtkCompositeTextureOnSurfaceFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);
}
