/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-18
 * \brief  Definition of the VtkCompositeLineToTubeFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkCompositeFilter.h"

/// @brief Converts lines to tube-objects.
class VtkCompositeLineToTubeFilter : public VtkCompositeFilter
{
public:
    VtkCompositeLineToTubeFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeLineToTubeFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

private:

};
