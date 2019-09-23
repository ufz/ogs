/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-19
 * \brief  Definition of the VtkCompositeImageToCylindersFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkCompositeFilter.h"

class VtkImageDataToLinePolyDataFilter;

/// @brief Creates cylinders that stand on top of the image with the length
/// of the corresponding first sub-pixel value (the red value). Useful to
/// visualize precipitation maps as a 3d bar chart.
class VtkCompositeImageToCylindersFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeImageToCylindersFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeImageToCylindersFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;

private:
    VtkImageDataToLinePolyDataFilter* _lineFilter;
};
