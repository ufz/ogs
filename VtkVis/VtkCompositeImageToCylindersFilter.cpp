/**
 * \file VtkCompositeImageToCylindersFilter.cpp
 * 19/10/2010 LB Initial implementation
 *
 * Implementation of VtkCompositeImageToCylindersFilter class
 */

// ** INCLUDES **
#include "VtkCompositeImageToCylindersFilter.h"

#include "VtkImageDataToLinePolyDataFilter.h"

#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>
#include <vtkUnsignedCharArray.h>

#include <QMap>
#include <QString>
#include <QVariant>
#include <QVector>

VtkCompositeImageToCylindersFilter::VtkCompositeImageToCylindersFilter(
        vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm)
{
	this->init();
}

void VtkCompositeImageToCylindersFilter::init()
{
	this->_inputDataObjectType = VTK_IMAGE_DATA;
	this->_outputDataObjectType = VTK_POLY_DATA;

	_lineFilter = VtkImageDataToLinePolyDataFilter::New();
	_lineFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
	_lineFilter->SetLengthScaleFactor(1);
	(*_algorithmUserProperties)["LengthScaleFactor"] = 1.0;
	_lineFilter->Update();

	vtkLookupTable* colormap = vtkLookupTable::New();
	colormap->SetTableRange(0, 100);
	colormap->SetHueRange(0.0, 0.666);
	colormap->SetNumberOfTableValues(256);
	colormap->ForceBuild();
	QList<QVariant> tableRangeList;
	tableRangeList.push_back(0);
	tableRangeList.push_back(100);
	QList<QVariant> hueRangeList;
	hueRangeList.push_back(0.0);
	hueRangeList.push_back(0.666);
	(*_algorithmUserVectorProperties)["TableRange"] = tableRangeList;
	(*_algorithmUserVectorProperties)["HueRange"] = hueRangeList;

	this->SetLookUpTable("Colors", colormap);

	vtkTubeFilter* tubeFilter = vtkTubeFilter::New();
	tubeFilter->SetInputConnection(_lineFilter->GetOutputPort());
	tubeFilter->CappingOn();
	tubeFilter->SetNumberOfSides(6);
	tubeFilter->SetRadius(_lineFilter->GetImageSpacing() * 0.25);
	(*_algorithmUserProperties)["NumberOfColors"] = 256;
	(*_algorithmUserProperties)["Capping"] = true;
	(*_algorithmUserProperties)["NumberOfSides"] = 6;
	(*_algorithmUserProperties)["RadiusFactor"] = 0.25;

	_outputAlgorithm = tubeFilter;
}

void VtkCompositeImageToCylindersFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	_lineFilter->SetUserProperty(name, value);

	// VtkImageDataToLinePolyDataFilter is equal to _firstAlgorithm
	// vtkTubeFilter is equal _outputAlgorithm
	if (name.compare("NumberOfColors") == 0)
	{
		vtkLookupTable* lut = this->GetLookupTable("Colors");
		if(lut)
			lut->SetNumberOfTableValues(value.toInt());
	}
	else if (name.compare("NumberOfSides") == 0)
		static_cast<vtkTubeFilter*>(_outputAlgorithm)->SetNumberOfSides(value.toInt());
	else if (name.compare("Capping") == 0)
		static_cast<vtkTubeFilter*>(_outputAlgorithm)->SetCapping(value.toBool());
	else if (name.compare("RadiusFactor") == 0)
		static_cast<vtkTubeFilter*>(_outputAlgorithm)->SetRadius(
		        _lineFilter->GetImageSpacing() * value.toDouble());
}

void VtkCompositeImageToCylindersFilter::SetUserVectorProperty( QString name,
                                                                QList<QVariant> values )
{
	VtkAlgorithmProperties::SetUserVectorProperty(name, values);

	_lineFilter->SetUserVectorProperty(name, values);

	if (name.compare("TableRange") == 0)
	{
		vtkLookupTable* lut = this->GetLookupTable("Colors");
		if(lut)
			lut->SetTableRange(values[0].toInt(), values[1].toInt());
	}
	else if (name.compare("HueRange") == 0)
	{
		vtkLookupTable* lut = this->GetLookupTable("Colors");
		if(lut)
			lut->SetHueRange(values[0].toDouble(), values[1].toDouble());
	}
}

VtkCompositeImageToCylindersFilter::~VtkCompositeImageToCylindersFilter()
{
	_lineFilter->Delete();
}
