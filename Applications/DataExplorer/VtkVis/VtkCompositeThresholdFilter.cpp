/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-25
 * \brief  Implementation of the VtkCompositeThresholdFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeThresholdFilter.h"

#include <logog/include/logog.hpp>

#include <vtkCellData.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>

VtkCompositeThresholdFilter::VtkCompositeThresholdFilter( vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm)
{
    this->init();
}

VtkCompositeThresholdFilter::~VtkCompositeThresholdFilter()
{
}

void VtkCompositeThresholdFilter::init()
{
    // Set meta information about input and output data types
    this->_inputDataObjectType = VTK_DATA_SET;
    this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

    // Because this is the only filter here we cannot use vtkSmartPointer
    vtkThreshold* threshold = vtkThreshold::New();
    threshold->SetInputConnection(_inputAlgorithm->GetOutputPort());

    // Use first array of parent as input array
    _inputAlgorithm->Update();
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(
        _inputAlgorithm->GetOutputDataObject(0));
    vtkDataSetAttributes* pointAttributes =
        dataSet->GetAttributes(vtkDataObject::AttributeTypes::POINT);
    vtkDataSetAttributes* cellAttributes =
        dataSet->GetAttributes(vtkDataObject::AttributeTypes::CELL);
    if(pointAttributes->GetNumberOfArrays() > 0)
        threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
            pointAttributes->GetArray(0)->GetName());
    else if(cellAttributes->GetNumberOfArrays() > 0)
        threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
            cellAttributes->GetArray(0)->GetName());
    else {
        WARN("Threshold filter could not find an array on its input object!")
        return;
    }

    // Sets a filter property which will be user editable
    threshold->SetSelectedComponent(0);

    // Setting the threshold to min / max values to ensure that the whole data
    // is first processed. This is needed for correct lookup table generation.
    const double dMin = std::numeric_limits<double>::lowest();
    const double dMax = std::numeric_limits<double>::max();
    threshold->ThresholdBetween(dMin, dMax);

    // Create a list for the ThresholdBetween (vector) property.
    QList<QVariant> thresholdRangeList;
    // Insert the values (same values as above)
    thresholdRangeList.push_back(dMin);
    thresholdRangeList.push_back(dMax);
    // Put that list in the property map
    (*_algorithmUserVectorProperties)["Range"] = thresholdRangeList;

    // Make a new entry in the property map for the SelectedComponent property
    (*_algorithmUserProperties)["Selected Component"] = 0;

    // Must all scalars match the criterium
    threshold->SetAllScalars(1);
    (*_algorithmUserProperties)["Evaluate all points"] = true;

    // The threshold filter is last one and so it is also the _outputAlgorithm
    _outputAlgorithm = threshold;
}

void VtkCompositeThresholdFilter::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);

    // Use the same name as in init()
    if (name.compare("Selected Component") == 0)
        // Set the property on the algorithm
        static_cast<vtkThreshold*>(_outputAlgorithm)->SetSelectedComponent(value.toInt());
    else if (name.compare("Evaluate all points") == 0)
        static_cast<vtkThreshold*>(_outputAlgorithm)->SetAllScalars(value.toBool());
}

void VtkCompositeThresholdFilter::SetUserVectorProperty( QString name, QList<QVariant> values )
{
    VtkAlgorithmProperties::SetUserVectorProperty(name, values);

    // Use the same name as in init()
    if (name.compare("Range") == 0)
        // Set the vector property on the algorithm
        static_cast<vtkThreshold*>(_outputAlgorithm)->ThresholdBetween(
                values[0].toDouble(), values[1].toDouble());
}
