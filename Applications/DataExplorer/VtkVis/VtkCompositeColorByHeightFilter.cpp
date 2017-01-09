/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Implementation of the VtkCompositeColorByHeightFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

    DataHolderLib::Color a{{0, 0, 255, 255}};    // blue
    DataHolderLib::Color b{{0, 255, 0, 255}};    // green
    DataHolderLib::Color c{{255, 255, 0, 255}};  // yellow
    DataHolderLib::Color d{{155, 100, 50, 255}}; // brown
    DataHolderLib::Color e{{255, 0, 0, 255}};    // red
    VtkColorLookupTable* ColorLookupTable = heightFilter->GetColorLookupTable();
    ColorLookupTable->setInterpolationType(DataHolderLib::LUTType::LINEAR);
    ColorLookupTable->setColor(-50, a);
    ColorLookupTable->setColor(0, a);
    ColorLookupTable->setColor(1, b);   // instant change at 0m a.s.l.
    ColorLookupTable->setColor(200, b); // green at about 200m a.s.l.
    ColorLookupTable->setColor(500, c); // yellow at about 500m and changing to red from then on
    ColorLookupTable->setColor(1000, d);
    ColorLookupTable->setColor(2000, e);
    ColorLookupTable->SetTableRange(-35, 2000);
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
