/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-05
 * \brief  Implementation of the VtkCompositeContourFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeContourFilter.h"

#include <vtkPointData.h>
#include <vtkContourFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <limits>

VtkCompositeContourFilter::VtkCompositeContourFilter( vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm)
{
    this->init();
}

VtkCompositeContourFilter::~VtkCompositeContourFilter() = default;

void VtkCompositeContourFilter::init()
{
    // Set meta information about input and output data types
    this->_inputDataObjectType = VTK_UNSTRUCTURED_GRID; //VTK_DATA_SET;
    this->_outputDataObjectType = VTK_UNSTRUCTURED_GRID;

    // Because this is the only filter here we cannot use vtkSmartPointer
    vtkContourFilter* contour = vtkContourFilter::New();
    contour->SetInputConnection(_inputAlgorithm->GetOutputPort());

    // Getting the scalar range from the active point data scalar of the input algorithm
    // This assumes that we do not want to contour on cell data.
    double range[2];
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(_inputAlgorithm->GetOutputDataObject(0));
    if(dataSet)
    {
        vtkPointData* pointData = dataSet->GetPointData();
        if(pointData)
            pointData->GetScalars()->GetRange(range);
    }
    else
    {
        // Setting the range to min / max values, this will result in a "bad table range"
        // vtk warning.
        range[0] = std::numeric_limits<double>::lowest();
        range[1] = std::numeric_limits<double>::max();
    }

    // Sets a filter vector property which will be user editable
    contour->GenerateValues(10, range[0], range[1]);

    // Create a list for the ThresholdBetween (vector) property.
    QList<QVariant> contourRangeList;
    // Insert the values (same values as above)
    contourRangeList.push_back(range[0]);
    contourRangeList.push_back(range[1]);
    // Put that list in the property map
    (*_algorithmUserVectorProperties)["Range"] = contourRangeList;

    // Make a new entry in the property map for the "Number of Values" property
    (*_algorithmUserProperties)["Number of Contours"] = 10;

    // The threshold filter is last one and so it is also the _outputAlgorithm
    _outputAlgorithm = contour;
}

void VtkCompositeContourFilter::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);

    // Use the same name as in init()
    if (name.compare("Number of Contours") == 0)
        static_cast<vtkContourFilter*>(_outputAlgorithm)->SetNumberOfContours(value.toInt());
}

void VtkCompositeContourFilter::SetUserVectorProperty( QString name, QList<QVariant> values )
{
    VtkAlgorithmProperties::SetUserVectorProperty(name, values);

    // Use the same name as in init()
    if (name.compare("Range") == 0)
        static_cast<vtkContourFilter*>(_outputAlgorithm)->GenerateValues(
                VtkAlgorithmProperties::GetUserProperty("Number of Contours").toInt(),
                values[0].toDouble(),
                values[1].toDouble());
}
