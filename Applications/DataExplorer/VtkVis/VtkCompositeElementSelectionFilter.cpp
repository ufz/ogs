/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-10
 * \brief  Implementation of the VtkCompositeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
: VtkCompositeFilter(inputAlgorithm), range_(0.0, 1.0), selection_name_("Selection")
{
}

void VtkCompositeElementSelectionFilter::init()
{
    double thresholdLower(range_.first);
    double thresholdUpper(range_.second);
    this->inputDataObjectType_ = VTK_UNSTRUCTURED_GRID;
    this->outputDataObjectType_ = VTK_UNSTRUCTURED_GRID;

    this->SetLookUpTable(QString::fromStdString(selection_name_), this->GetLookupTable());
    vtkSmartPointer<VtkAppendArrayFilter> selFilter(nullptr);
    if (!selection_.empty())
    {
        selFilter = vtkSmartPointer<VtkAppendArrayFilter>::New();
        selFilter->SetInputConnection(inputAlgorithm_->GetOutputPort());
        selFilter->SetArray(selection_name_, selection_);
        selFilter->Update();
    }

    vtkSmartPointer<vtkIdFilter> idFilter = vtkSmartPointer<vtkIdFilter>::New();
    if (selection_.empty())
    {  // if the array is empty it is assumed that an existing array should be
       // used
        idFilter->SetInputConnection(inputAlgorithm_->GetOutputPort());
    }
    else
    {
        idFilter->SetInputConnection(selFilter->GetOutputPort());
    }
        idFilter->PointIdsOn();
        idFilter->CellIdsOn();
        idFilter->FieldDataOn();
        idFilter->Update();

    vtkThreshold* threshold = vtkThreshold::New();
        threshold->SetInputConnection(idFilter->GetOutputPort());
        threshold->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS, selection_name_.c_str());
        threshold->SetSelectedComponent(0);
        threshold->ThresholdBetween(thresholdLower, thresholdUpper);
        threshold->Update();

    QList<QVariant> thresholdRangeList;
    thresholdRangeList.push_back(thresholdLower);
    thresholdRangeList.push_back(thresholdUpper);
    (*algorithmUserVectorProperties_)["Threshold Between"] = thresholdRangeList;
    outputAlgorithm_ = threshold;
}

void VtkCompositeElementSelectionFilter::setSelectionArray(const std::string &selection_name, const std::vector<double> &selection)
{
    selection_name_ = selection_name;
    selection_ = selection;
    init();
}

void VtkCompositeElementSelectionFilter::SetUserVectorProperty( QString name, QList<QVariant> values)
{
    VtkAlgorithmProperties::SetUserVectorProperty(name, values);

    if (name.compare("Threshold Between") == 0)
    {
        static_cast<vtkThreshold*>(outputAlgorithm_)
            ->ThresholdBetween(values[0].toDouble(), values[1].toDouble());
    }
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
