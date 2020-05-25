/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkCompositeImageToPointCloudFilter.h"

#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <QMap>
#include <QString>
#include <QVariant>

#include "VtkImageDataToPointCloudFilter.h"

VtkCompositeImageToPointCloudFilter::VtkCompositeImageToPointCloudFilter(vtkAlgorithm* inputAlgorithm)
    : VtkCompositeFilter(inputAlgorithm)
{
    init();
}

void VtkCompositeImageToPointCloudFilter::init()
{
    inputDataObjectType_ = VTK_IMAGE_DATA;
    outputDataObjectType_ = VTK_POLY_DATA;

    VtkImageDataToPointCloudFilter* point_cloud_filter =
        VtkImageDataToPointCloudFilter::New();
    point_cloud_filter->SetInputConnection(inputAlgorithm_->GetOutputPort());
    inputAlgorithm_->Update();

    QList<QVariant> n_points_range_list;
    n_points_range_list.push_back(point_cloud_filter->GetMinNumberOfPointsPerCell());
    n_points_range_list.push_back(point_cloud_filter->GetMaxNumberOfPointsPerCell());
    (*algorithmUserVectorProperties_)["Number of points range"] = n_points_range_list;
    QList<QVariant> vertical_extent_list;
    vertical_extent_list.push_back(point_cloud_filter->GetMinHeight());
    vertical_extent_list.push_back(point_cloud_filter->GetMaxHeight());
    (*algorithmUserVectorProperties_)["Vertical extent"] = vertical_extent_list;
    (*algorithmUserProperties_)["Logarithmic interpolation"] = !point_cloud_filter->GetIsLinear();
    (*algorithmUserProperties_)["Gamma value"] = point_cloud_filter->GetGamma();

    point_cloud_filter->Update();
    outputAlgorithm_ = point_cloud_filter;
}

void VtkCompositeImageToPointCloudFilter::SetUserProperty(QString name, QVariant value)
{
    VtkAlgorithmProperties::SetUserProperty(name, value);

    if ((name == "Gamma value") && (value.toDouble() > 0))
    {
        static_cast<VtkImageDataToPointCloudFilter*>(outputAlgorithm_)->SetGamma(value.toDouble());
    }
    if (name == "Logarithmic interpolation")
    {
        if (value.toBool())
        {
            double const gamma =
                VtkAlgorithmProperties::GetUserProperty("Gamma value").toDouble();
            if (gamma > 0)
            {
                static_cast<VtkImageDataToPointCloudFilter*>(outputAlgorithm_)
                    ->useLogarithmicInterpolation(gamma);
            }
        }
        else
        {
            static_cast<VtkImageDataToPointCloudFilter*>(outputAlgorithm_)
                ->useLinearInterpolation();
        }
    }
}

void VtkCompositeImageToPointCloudFilter::SetUserVectorProperty(QString name, QList<QVariant> values)
{
    VtkAlgorithmProperties::SetUserVectorProperty(name, values);

    if (name == "Number of points range")
    {
        if (values[0].toInt() >= 0 && values[1].toInt() >= 0 &&
            values[0].toInt() <= values[1].toInt())
        {
            static_cast<VtkImageDataToPointCloudFilter*>(outputAlgorithm_)
                ->SetMinNumberOfPointsPerCell(values[0].toInt());
            static_cast<VtkImageDataToPointCloudFilter*>(outputAlgorithm_)
                ->SetMaxNumberOfPointsPerCell(values[1].toInt());
        }
    }
    else if (name == "Vertical extent")
    {
        if (values[0].toDouble() <= values[1].toDouble())
        {
            static_cast<VtkImageDataToPointCloudFilter*>(outputAlgorithm_)
                ->SetMinHeight(values[0].toDouble());
            static_cast<VtkImageDataToPointCloudFilter*>(outputAlgorithm_)
                ->SetMaxHeight(values[1].toDouble());
        }
    }
}
