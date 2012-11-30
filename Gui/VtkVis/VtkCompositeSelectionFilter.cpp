/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkCompositeSelectionFilter.cpp
 *
 * Created on 2011-02-10 by Karsten Rink
 */

// ** INCLUDES **
#include "VtkCompositeSelectionFilter.h"
#include "VtkSelectionFilter.h"
#include "VtkColorLookupTable.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkIdFilter.h>
#include <vtkUnstructuredGrid.h>

VtkCompositeSelectionFilter::VtkCompositeSelectionFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm)
{
	//this->init();
}

void VtkCompositeSelectionFilter::init()
{
	const char* filter_name("Selection");
	double thresholdLower(0.0), thresholdUpper(1.0);
	this->_inputDataObjectType = VTK_UNSTRUCTURED_GRID;
	this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

	this->SetLookUpTable(QString(filter_name), this->GetLookupTable());

	VtkSelectionFilter* selFilter = VtkSelectionFilter::New();
	selFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
	selFilter->SetSelectionArray(_selection, thresholdLower, thresholdUpper);
	selFilter->Update();

	vtkIdFilter* idFilter = vtkIdFilter::New();
		idFilter->SetInputConnection(selFilter->GetOutputPort());
		idFilter->PointIdsOn();
		idFilter->CellIdsOn();
		idFilter->FieldDataOn();
		idFilter->Update();

	vtkThreshold* threshold = vtkThreshold::New();
	threshold->SetInputConnection(idFilter->GetOutputPort());
	threshold->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS, filter_name);
	threshold->SetSelectedComponent(0);
	threshold->ThresholdBetween(thresholdLower, thresholdUpper);
	threshold->Update();

	QList<QVariant> thresholdRangeList;
	thresholdRangeList.push_back(0.0);
	thresholdRangeList.push_back(1.0);
	(*_algorithmUserVectorProperties)["Threshold Between"] = thresholdRangeList;

	_outputAlgorithm = threshold;
}

void VtkCompositeSelectionFilter::SetUserVectorProperty( QString name, QList<QVariant> values)
{
	VtkAlgorithmProperties::SetUserVectorProperty(name, values);

	if (name.compare("Threshold Between") == 0)
		static_cast<vtkThreshold*>(_outputAlgorithm)->ThresholdBetween(
		        values[0].toDouble(), values[1].toDouble());
}

VtkColorLookupTable* VtkCompositeSelectionFilter::GetLookupTable()
{
	VtkColorLookupTable* lut = VtkColorLookupTable::New();
	lut->SetTableRange(0,1);
	unsigned char a[4] = { 0, 0, 255, 255 }; // blue
	unsigned char b[4] = { 0, 255, 0, 255 }; // green
	unsigned char c[4] = { 255, 255, 0, 255 }; // yellow
	unsigned char d[4] = { 255, 0, 0, 255 }; // red
	lut->setColor(1.0, a);
	lut->setColor(0.5, b);
	lut->setColor(0.25, c);
	lut->setColor(0.1, d);
	lut->Build();
	return lut;
}
