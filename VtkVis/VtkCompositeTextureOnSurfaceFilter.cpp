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

#include "OGSRaster.h"

#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>

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
	                                                settings.value(
	                                                        "lastOpenedTextureFileDirectory").
	                                                toString(),
	                                                "Raster files (*.asc *.bmp *.jpg *.png *.tif);;");
	QFileInfo fi(fileName);

	if ((fi.suffix().toLower() == "asc") || (fi.suffix().toLower() == "tif") ||
	    (fi.suffix().toLower() == "png") ||
	    (fi.suffix().toLower() == "jpg") || (fi.suffix().toLower() == "bmp"))
	{
		QImage img;
		QPointF origin;
		double scalingFactor = 0;

		OGSRaster::loadImage(fileName, img, origin, scalingFactor);
		std::pair<float, float> org(origin.x(), origin.y());
		surface->SetRaster(img, org, scalingFactor);
		surface->Update();

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedTextureFileDirectory", dir.absolutePath());
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
