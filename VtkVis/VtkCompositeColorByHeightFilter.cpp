/**
 * \file VtkCompositeColorByHeightFilter.cpp
 * 01/11/2010 KR Initial implementation
 *
 * Implementation of VtkCompositePointToGlyphFilter class
 */

// ** INCLUDES **
#include "VtkColorByHeightFilter.h"
#include "VtkColorLookupTable.h"
#include "VtkCompositeColorByHeightFilter.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

VtkCompositeColorByHeightFilter::VtkCompositeColorByHeightFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm)
{
	this->init();
}

void VtkCompositeColorByHeightFilter::init()
{
	this->_inputDataObjectType = VTK_DATA_SET;
	this->_outputDataObjectType = VTK_POLY_DATA;

	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter;
	VtkColorByHeightFilter* heightFilter = VtkColorByHeightFilter::New();

	if (dynamic_cast<vtkUnstructuredGrid*>(_inputAlgorithm->GetOutputDataObject(0)))
	{
		surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
		surfaceFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
		heightFilter->SetInputConnection(surfaceFilter->GetOutputPort());
	}
	else
		heightFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());

	heightFilter->SetTableRange(-35, 800); // default min- and max-height (default: blue to red)
	//heightFilter->GetColorLookupTable()->setInterpolationType(ColorLookupTable::EXPONENTIAL);
	unsigned char a[4] = { 0, 0, 255, 255 }; // blue
	unsigned char b[4] = { 0, 255, 0, 255 }; // green
	unsigned char c[4] = { 255, 255, 0, 255 }; // yellow
	unsigned char d[4] = { 255, 0, 0, 255 }; // red
	//unsigned char e[4] = { 255, 255, 255, 255 }; // white
	heightFilter->GetColorLookupTable()->setColor(0.0, a);
	heightFilter->GetColorLookupTable()->setColor(0.2, b); // green at about 150m
	heightFilter->GetColorLookupTable()->setColor(0.6, c); // yellow at about 450m and changing to red from then on
	heightFilter->GetColorLookupTable()->setColor(1.0, d);
	//heightFilter->GetColorLookupTable()->setColor(1.0, e);
	heightFilter->Update();

	_outputAlgorithm = heightFilter;
	_activeAttributeName = heightFilter->GetActiveAttribute();
}

void VtkCompositeColorByHeightFilter::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}
