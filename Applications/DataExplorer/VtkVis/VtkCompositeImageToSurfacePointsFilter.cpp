/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkCompositeImageToSurfacePointsFilter.h"

#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <QMap>
#include <QString>
#include <QVariant>

#include "VtkImageDataToSurfacePointsFilter.h"

VtkCompositeImageToSurfacePointsFilter::VtkCompositeImageToSurfacePointsFilter(
    vtkAlgorithm* inputAlgorithm)
    : VtkCompositeFilter(inputAlgorithm)
{
    init();
}

void VtkCompositeImageToSurfacePointsFilter::init()
{
    inputDataObjectType_ = VTK_IMAGE_DATA;
    outputDataObjectType_ = VTK_POLY_DATA;

    VtkImageDataToSurfacePointsFilter* point_cloud_filter =
        VtkImageDataToSurfacePointsFilter::New();
    point_cloud_filter->SetInputConnection(inputAlgorithm_->GetOutputPort());
    inputAlgorithm_->Update();

    (*algorithmUserProperties_)["Points per pixel"] =
        static_cast<int>(point_cloud_filter->GetPointsPerPixel());
    point_cloud_filter->Update();
    outputAlgorithm_ = point_cloud_filter;
}

void VtkCompositeImageToSurfacePointsFilter::SetUserProperty(QString name, QVariant value)
{
    VtkAlgorithmProperties::SetUserProperty(name, value);
    if ((name == "Points per pixel") && (value.toInt() > 0))
    {
        static_cast<VtkImageDataToSurfacePointsFilter*>(outputAlgorithm_)
            ->SetPointsPerPixel(static_cast<vtkIdType>(value.toInt()));
    }
}

void VtkCompositeImageToSurfacePointsFilter::SetUserVectorProperty(
    QString name, QList<QVariant> values)
{
    Q_UNUSED(name);
    Q_UNUSED(values);
}
