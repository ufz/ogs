/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-10
 * \brief  Implementation of the VtkCompositeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeElementSelectionFilter.h"
#include "VtkAppendArrayFilter.h"
#include "VtkCompositePointToGlyphFilter.h"
#include "VtkColorLookupTable.h"
#include "VtkPointsSource.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkIdFilter.h>
#include <vtkUnstructuredGrid.h>

#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkPointData.h>


VtkCompositeElementSelectionFilter::VtkCompositeElementSelectionFilter( vtkAlgorithm* inputAlgorithm )
: VtkCompositeFilter(inputAlgorithm), _range(0.0, 1.0), _selection_name("Selection")
{
}

void VtkCompositeElementSelectionFilter::init()
{
    double thresholdLower(_range.first), thresholdUpper(_range.second);
    this->_inputDataObjectType = VTK_UNSTRUCTURED_GRID;
    this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

    this->SetLookUpTable(QString::fromStdString(_selection_name), this->GetLookupTable());
    vtkSmartPointer<VtkAppendArrayFilter> selFilter(nullptr);
    if (!_selection.empty())
    {
        selFilter = vtkSmartPointer<VtkAppendArrayFilter>::New();
        selFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
        selFilter->SetArray(_selection_name, _selection);
        selFilter->Update();
    }

    vtkSmartPointer<vtkIdFilter> idFilter = vtkSmartPointer<vtkIdFilter>::New();
        if (_selection.empty()) // if the array is empty it is assumed that an existing array should be used
            idFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
        else
            idFilter->SetInputConnection(selFilter->GetOutputPort());
        idFilter->PointIdsOn();
        idFilter->CellIdsOn();
        idFilter->FieldDataOn();
        idFilter->Update();

    vtkThreshold* threshold = vtkThreshold::New();
        threshold->SetInputConnection(idFilter->GetOutputPort());
        threshold->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS, _selection_name.c_str());
        threshold->SetSelectedComponent(0);
        threshold->ThresholdBetween(thresholdLower, thresholdUpper);
        threshold->Update();

    QList<QVariant> thresholdRangeList;
    thresholdRangeList.push_back(thresholdLower);
    thresholdRangeList.push_back(thresholdUpper);
    (*_algorithmUserVectorProperties)["Threshold Between"] = thresholdRangeList;
    _outputAlgorithm = threshold;
}

void VtkCompositeElementSelectionFilter::setSelectionArray(const std::string &selection_name, const std::vector<double> &selection)
{
    _selection_name = selection_name;
    _selection = selection;
    init();
}

void VtkCompositeElementSelectionFilter::SetUserVectorProperty( QString name, QList<QVariant> values)
{
    VtkAlgorithmProperties::SetUserVectorProperty(name, values);

    if (name.compare("Threshold Between") == 0)
        static_cast<vtkThreshold*>(_outputAlgorithm)->ThresholdBetween(
                values[0].toDouble(), values[1].toDouble());
}

VtkColorLookupTable* VtkCompositeElementSelectionFilter::GetLookupTable()
{
    VtkColorLookupTable* lut = VtkColorLookupTable::New();
    lut->SetTableRange(0,1);
    DataHolderLib::Color a{{0, 0, 255, 255}};   // blue
    DataHolderLib::Color b{{0, 255, 0, 255}};   // green
    DataHolderLib::Color c{{255, 255, 0, 255}};  // yellow
    DataHolderLib::Color d{{255, 0, 0, 255}};    // red
    lut->setColor(1.0, a);
    lut->setColor(0.5, b);
    lut->setColor(0.25, c);
    lut->setColor(0.1, d);
    lut->Build();
    return lut;
}
