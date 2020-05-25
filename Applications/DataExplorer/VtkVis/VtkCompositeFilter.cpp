/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-19
 * \brief  Implementation of the VtkCompositeFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeFilter.h"

#include <vtkAlgorithm.h>
#include <vtkPolyData.h>

#include <QMapIterator>
#include <QString>
#include <QVector>

VtkCompositeFilter::VtkCompositeFilter(vtkAlgorithm* inputAlgorithm)
    : inputDataObjectType_(0), outputDataObjectType_(1),
      inputAlgorithm_(inputAlgorithm), outputAlgorithm_(nullptr)
{
}

VtkCompositeFilter::~VtkCompositeFilter()
{
    outputAlgorithm_->Delete();
}

double VtkCompositeFilter::GetInitialRadius() const
{
    double bounding_box[6];
    static_cast<vtkPolyData*>(this->inputAlgorithm_->GetOutputDataObject(0))->GetBounds(bounding_box);
    double x_diff = fabs(bounding_box[0]-bounding_box[1]);
    double y_diff = fabs(bounding_box[2]-bounding_box[3]);
    double z_diff = fabs(bounding_box[4]-bounding_box[5]);

    double max = (x_diff == 0) ? 1 : x_diff;
    max = (max > y_diff) ? max : y_diff;
    max = (max > z_diff) ? max : z_diff;

    return max/200.0;
}
