/**
 * \file VtkCompositeColormapToImageFilter.cpp
 * 21/10/2010 LB Initial implementation
 * 
 * Implementation of VtkCompositeColormapToImageFilter class
 */

// ** INCLUDES **
#include "VtkCompositeColormapToImageFilter.h"

#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>
#include <vtkSmartPointer.h>

VtkCompositeColormapToImageFilter::VtkCompositeColormapToImageFilter( vtkAlgorithm* inputAlgorithm )
: VtkCompositeFilter(inputAlgorithm)
{
	this->init();
}

VtkCompositeColormapToImageFilter::~VtkCompositeColormapToImageFilter()
{
}

void VtkCompositeColormapToImageFilter::init()
{
	this->_inputDataObjectType = VTK_IMAGE_DATA;
	this->_outputDataObjectType = VTK_IMAGE_DATA;

	vtkSmartPointer<vtkLookupTable> colormap = vtkSmartPointer<vtkLookupTable>::New();
	colormap->SetTableRange(0, 100);
	colormap->SetHueRange(0.0, 0.666);
	colormap->SetNumberOfTableValues(256);
	QList<QVariant> tableRangeList;
	tableRangeList.push_back(0);
	tableRangeList.push_back(100);
	QList<QVariant> hueRangeList;
	hueRangeList.push_back(0.0);
	hueRangeList.push_back(0.666);
	(*_algorithmUserVectorProperties)["TableRange"] = tableRangeList;
	(*_algorithmUserVectorProperties)["HueRange"] = hueRangeList;

	vtkImageMapToColors* map = vtkImageMapToColors::New();
	map->SetInputConnection(0, _inputAlgorithm->GetOutputPort());
	map->SetLookupTable(colormap);
	map->SetPassAlphaToOutput(1);
	(*_algorithmUserProperties)["PassAlphaToOutput"] = true;
	(*_algorithmUserProperties)["NumberOfColors"] = 256;

	_outputAlgorithm = map;
}

void VtkCompositeColormapToImageFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	vtkImageMapToColors* map = static_cast<vtkImageMapToColors*>(_outputAlgorithm);
	if (name.compare("PassAlphaToOutput") == 0)
		map->SetPassAlphaToOutput(value.toBool());
	else if (name.compare("NumberOfColors") == 0)
		static_cast<vtkLookupTable*>(map->GetLookupTable())->SetNumberOfTableValues(value.toInt());
}

void VtkCompositeColormapToImageFilter::SetUserVectorProperty( QString name, QList<QVariant> values )
{
	VtkAlgorithmProperties::SetUserVectorProperty(name, values);

	vtkImageMapToColors* map = static_cast<vtkImageMapToColors*>(_outputAlgorithm);
	if (name.compare("TableRange") == 0)
		static_cast<vtkLookupTable*>(map->GetLookupTable())->SetTableRange(values[0].toInt(), values[1].toInt());
	else if (name.compare("HueRange") == 0)
		static_cast<vtkLookupTable*>(map->GetLookupTable())->SetHueRange(values[0].toDouble(), values[1].toDouble());
}
