/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Implementation of the VtkCompositeColorByHeightFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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

	unsigned char a[4] = { 0, 0, 255, 255 }; // blue
	unsigned char b[4] = { 0, 255, 0, 255 }; // green
	unsigned char c[4] = { 255, 255, 0, 255 }; // yellow
	unsigned char d[4] = { 255, 0, 0, 255 }; // red
	VtkColorLookupTable* ColorLookupTable = heightFilter->GetColorLookupTable();
	ColorLookupTable->setInterpolationType(VtkColorLookupTable::LUTType::LINEAR);
	ColorLookupTable->setColor(-35, a);
	ColorLookupTable->setColor(150, b); // green at about 150m
	ColorLookupTable->setColor(450, c); // yellow at about 450m and changing to red from then on
	ColorLookupTable->setColor(800, d);
	ColorLookupTable->SetTableRange(-35, 800);
	ColorLookupTable->Build();

	// This passes ownership of the ColorLookupTable to VtkVisPointSetItem
	heightFilter->SetLookUpTable("P-Colors", ColorLookupTable);
	heightFilter->Update();

	_outputAlgorithm = heightFilter;
	_activeAttributeName = heightFilter->GetActiveAttribute();
}

void VtkCompositeColorByHeightFilter::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}
