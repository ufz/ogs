/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

VtkCompositeImageToSurfacePointsFilter::~VtkCompositeImageToSurfacePointsFilter() = default;

void VtkCompositeImageToSurfacePointsFilter::init()
{
    _inputDataObjectType = VTK_IMAGE_DATA;
    _outputDataObjectType = VTK_POLY_DATA;

    VtkImageDataToSurfacePointsFilter* point_cloud_filter =
        VtkImageDataToSurfacePointsFilter::New();
    point_cloud_filter->SetInputConnection(_inputAlgorithm->GetOutputPort());
    _inputAlgorithm->Update();

    (*_algorithmUserProperties)["Points per pixel"] =
        static_cast<int>(point_cloud_filter->GetPointsPerPixel());
    point_cloud_filter->Update();
    _outputAlgorithm = point_cloud_filter;
}

void VtkCompositeImageToSurfacePointsFilter::SetUserProperty(QString name, QVariant value)
{
    VtkAlgorithmProperties::SetUserProperty(name, value);
    if ((name == "Points per pixel") && (value.toInt() > 0))
        static_cast<VtkImageDataToSurfacePointsFilter*>(_outputAlgorithm)
            ->SetPointsPerPixel(static_cast<vtkIdType>(value.toInt()));
}

void VtkCompositeImageToSurfacePointsFilter::SetUserVectorProperty(
    QString name, QList<QVariant> values)
{
    Q_UNUSED(name);
    Q_UNUSED(values);
}
