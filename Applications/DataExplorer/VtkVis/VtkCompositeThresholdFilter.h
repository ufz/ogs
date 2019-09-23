/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-15
 * \brief  Definition of the VtkCompositeThresholdFilter class.
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

/// @brief Visualises only parts of meshes that are above/below/within given thresholds.
/// In init() the threshold is first set to double min / max values. Set the
/// threshold later on via SetUserVectorProperty() to the actual data range.
class VtkCompositeThresholdFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeThresholdFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeThresholdFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;

private:
};
