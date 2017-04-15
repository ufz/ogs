/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-05
 * \brief  Definition of the VtkCompositeContourFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkCompositeFilter.h"

/// @brief Visualisation of contour-lines/-planes within dense scalar fields.
/// In init() the threshold is first set to double min / max values. Set the
/// threshold later on via SetUserVectorProperty() to the actual data range.
class VtkCompositeContourFilter : public VtkCompositeFilter
{
public:
    VtkCompositeContourFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeContourFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;

private:
};
