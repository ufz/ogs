/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Definition of the VtkCompositeColorByHeightFilter class.
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

class vtkSphereSource;

/// @brief This filter colors the input by the points z-value.
class VtkCompositeColorByHeightFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeColorByHeightFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeColorByHeightFilter() override = default;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

protected:
};
