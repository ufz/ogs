/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-19
 * \brief  Implementation of the VtkCompositeImageToCylindersFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeImageToCylindersFilter.h"

#include "VtkImageDataToLinePolyDataFilter.h"

#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>
#include <vtkUnsignedCharArray.h>

#include <QMap>
#include <QString>
#include <QVariant>
#include <QVector>

VtkCompositeImageToCylindersFilter::VtkCompositeImageToCylindersFilter(
        vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm), lineFilter_(nullptr)
{
    this->init();
}

void VtkCompositeImageToCylindersFilter::init()
{
    this->inputDataObjectType_ = VTK_IMAGE_DATA;
    this->outputDataObjectType_ = VTK_POLY_DATA;

    lineFilter_ = VtkImageDataToLinePolyDataFilter::New();
    lineFilter_->SetInputConnection(inputAlgorithm_->GetOutputPort());
    lineFilter_->SetLengthScaleFactor(1);
    (*algorithmUserProperties_)["LengthScaleFactor"] = 1.0;
    lineFilter_->Update();

    double range[2];
    // The data is always on points
    vtkDataSet::SafeDownCast(lineFilter_->GetOutputDataObject(0))->GetPointData()->GetScalars()->GetRange(range);

    vtkLookupTable* colormap = vtkLookupTable::New();
    colormap->SetTableRange(range[0], range[1]);
    colormap->SetHueRange(0.0, 0.666);
    colormap->SetNumberOfTableValues(256);
    colormap->ForceBuild();
    QList<QVariant> tableRangeList;
    tableRangeList.push_back(range[0]);
    tableRangeList.push_back(range[1]);
    QList<QVariant> hueRangeList;
    hueRangeList.push_back(0.0);
    hueRangeList.push_back(0.666);
    (*algorithmUserVectorProperties_)["TableRange"] = tableRangeList;
    (*algorithmUserVectorProperties_)["HueRange"] = hueRangeList;

    this->SetLookUpTable("P-Colors", colormap);

    vtkTubeFilter* tubeFilter = vtkTubeFilter::New();
    tubeFilter->SetInputConnection(lineFilter_->GetOutputPort());
    tubeFilter->CappingOn();
    tubeFilter->SetNumberOfSides(6);
    tubeFilter->SetRadius(lineFilter_->GetImageSpacing() * 0.25);
    (*algorithmUserProperties_)["NumberOfColors"] = 256;
    (*algorithmUserProperties_)["Capping"] = true;
    (*algorithmUserProperties_)["NumberOfSides"] = 6;
    (*algorithmUserProperties_)["RadiusFactor"] = 0.25;

    outputAlgorithm_ = tubeFilter;
}

void VtkCompositeImageToCylindersFilter::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);

    lineFilter_->SetUserProperty(name, value);

    // VtkImageDataToLinePolyDataFilter is equal to firstAlgorithm_
    // vtkTubeFilter is equal outputAlgorithm_
    if (name.compare("NumberOfColors") == 0)
    {
        vtkLookupTable* lut = this->GetLookupTable("P-Colors");
        if (lut)
        {
            lut->SetNumberOfTableValues(value.toInt());
        }
    }
    else if (name.compare("NumberOfSides") == 0)
    {
        static_cast<vtkTubeFilter*>(outputAlgorithm_)->SetNumberOfSides(value.toInt());
    }
    else if (name.compare("Capping") == 0)
    {
        static_cast<vtkTubeFilter*>(outputAlgorithm_)->SetCapping(value.toBool());
    }
    else if (name.compare("RadiusFactor") == 0)
    {
        static_cast<vtkTubeFilter*>(outputAlgorithm_)
            ->SetRadius(lineFilter_->GetImageSpacing() * value.toDouble());
    }
}

void VtkCompositeImageToCylindersFilter::SetUserVectorProperty( QString name,
                                                                QList<QVariant> values )
{
    VtkAlgorithmProperties::SetUserVectorProperty(name, values);

    lineFilter_->SetUserVectorProperty(name, values);

    if (name.compare("TableRange") == 0)
    {
        vtkLookupTable* lut = this->GetLookupTable("P-Colors");
        if (lut)
        {
            lut->SetTableRange(values[0].toDouble(), values[1].toDouble());
        }
    }
    else if (name.compare("HueRange") == 0)
    {
        vtkLookupTable* lut = this->GetLookupTable("P-Colors");
        if (lut)
        {
            lut->SetHueRange(values[0].toDouble(), values[1].toDouble());
        }
    }
}

VtkCompositeImageToCylindersFilter::~VtkCompositeImageToCylindersFilter()
{
    lineFilter_->Delete();
}
