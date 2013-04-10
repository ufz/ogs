/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-10
 * \brief  Implementation of the VtkCompositeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeSelectionFilter.h"
#include "VtkAppendArrayFilter.h"
#include "VtkCompositePointToGlyphFilter.h"
#include "VtkColorLookupTable.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkIdFilter.h>
#include <vtkUnstructuredGrid.h>

#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkPointData.h>


VtkCompositeSelectionFilter::VtkCompositeSelectionFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm), _selection_name("Selection"), _is_element_array(true)
{
	//this->init();
}

void VtkCompositeSelectionFilter::init()
{
	double thresholdLower(0.0), thresholdUpper(100000.0);
	this->_inputDataObjectType = VTK_UNSTRUCTURED_GRID;
	this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

	this->SetLookUpTable(QString::fromStdString(_selection_name), this->GetLookupTable());

	VtkAppendArrayFilter* selFilter (NULL);
	if (!_selection.empty())
	{
		selFilter = VtkAppendArrayFilter::New();
		selFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
		selFilter->SetArray(_selection_name, _selection);
		selFilter->Update();
	}

	vtkIdFilter* idFilter = vtkIdFilter::New();
		if (_selection.empty()) // if the array is empty it is assumed that an existing array should be used
			idFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
		else
			idFilter->SetInputConnection(selFilter->GetOutputPort());
		idFilter->PointIdsOn();
		idFilter->CellIdsOn();
		idFilter->FieldDataOn();
		idFilter->Update();
		
	size_t nPoints = idFilter->GetOutput()->GetPointData()->GetNumberOfTuples();

	_threshold = vtkThreshold::New();
		_threshold->SetInputConnection(idFilter->GetOutputPort());
		if (_is_element_array)
			_threshold->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS, _selection_name.c_str());
		else
			_threshold->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, _selection_name.c_str());
		_threshold->SetSelectedComponent(0);
		_threshold->ThresholdBetween(thresholdLower, thresholdUpper);
		_threshold->Update();
		
	QList<QVariant> thresholdRangeList;
	thresholdRangeList.push_back(thresholdLower);
	thresholdRangeList.push_back(thresholdUpper);
	(*_algorithmUserVectorProperties)["Threshold Between"] = thresholdRangeList;
	
	if (!_is_element_array)
	{
			size_t nPoints = _threshold->GetOutput()->GetPointData()->GetNumberOfTuples();
		VtkCompositeFilter* composite = new VtkCompositePointToGlyphFilter(_threshold);
		composite->SetUserProperty("Radius", this->GetInitialRadius());
		_outputAlgorithm = composite->GetOutputAlgorithm();
	}
	else
		_outputAlgorithm = _threshold;
}

void VtkCompositeSelectionFilter::setSelectionArray(const std::string &selection_name, bool is_element_array, const std::vector<double> &selection)
{
	_selection_name = selection_name;
	_is_element_array = is_element_array;
	_selection = selection;
	init(); 
}

void VtkCompositeSelectionFilter::SetUserVectorProperty( QString name, QList<QVariant> values)
{
	VtkAlgorithmProperties::SetUserVectorProperty(name, values);
	double a = values[0].toDouble();
	double b = values[1].toDouble();

	if (name.compare("Threshold Between") == 0)
		static_cast<vtkThreshold*>(_threshold)->ThresholdBetween(
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
