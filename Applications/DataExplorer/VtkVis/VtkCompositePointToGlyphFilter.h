/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-21
 * \brief  Definition of the VtkCompositePointToGlyphFilter class.
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

class vtkSphereSource;

/// @brief Converts point data to scalar-scaled spheres.
class VtkCompositePointToGlyphFilter : public VtkCompositeFilter
{
public:
    VtkCompositePointToGlyphFilter(vtkAlgorithm* inputAlgorithm);
    virtual ~VtkCompositePointToGlyphFilter();

    virtual void init() override;

    virtual void SetUserProperty(QString name, QVariant value) override;

private:
    vtkSphereSource* _glyphSource;
};
