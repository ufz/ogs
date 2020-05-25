/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-18
 * \brief  Implementation of the VtkCompositeLineToTubeFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeLineToTubeFilter.h"

#include <vtkCleanPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>

VtkCompositeLineToTubeFilter::VtkCompositeLineToTubeFilter( vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm)
{
    this->init();
}

VtkCompositeLineToTubeFilter::~VtkCompositeLineToTubeFilter() = default;

void VtkCompositeLineToTubeFilter::init()
{
    this->inputDataObjectType_ = VTK_DATA_SET;
    this->outputDataObjectType_ = VTK_POLY_DATA;

    // collapse coincident points
    vtkSmartPointer<vtkCleanPolyData> mergePoints = vtkSmartPointer<vtkCleanPolyData>::New();
    mergePoints->SetInputConnection(0, inputAlgorithm_->GetOutputPort(0));
    mergePoints->SetTolerance(0.0);
    mergePoints->ConvertLinesToPointsOn();

    double default_radius(GetInitialRadius());
    int default_number_of_sides(8);
    vtkTubeFilter* tubes = vtkTubeFilter::New();
    tubes->SetInputConnection(0, mergePoints->GetOutputPort(0));

    //tubes->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"StationValue");
    //tubes->SetVaryRadiusToVaryRadiusByScalar(); // KR radius changes with scalar

    tubes->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"Stratigraphies");
    tubes->SetRadius(default_radius);
    tubes->SetNumberOfSides(default_number_of_sides);
    tubes->SetCapping(1);

    (*algorithmUserProperties_)["Radius"] = default_radius;
    (*algorithmUserProperties_)["NumberOfSides"] = default_number_of_sides;
    (*algorithmUserProperties_)["Capping"] = true;

    outputAlgorithm_ = tubes;
}

void VtkCompositeLineToTubeFilter::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);

    if (name.compare("Radius") == 0)
    {
        static_cast<vtkTubeFilter*>(outputAlgorithm_)->SetRadius(value.toDouble());
    }
    else if (name.compare("NumberOfSides") == 0)
    {
        static_cast<vtkTubeFilter*>(outputAlgorithm_)->SetNumberOfSides(value.toInt());
    }
    else if (name.compare("Capping") == 0)
    {
        static_cast<vtkTubeFilter*>(outputAlgorithm_)
            ->SetCapping(value.toBool());
    }
}
