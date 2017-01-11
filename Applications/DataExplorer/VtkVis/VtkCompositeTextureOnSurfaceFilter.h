/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-18
 * \brief  Definition of the VtkCompositeTextureOnSurfaceFilter class.
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

/// @brief Puts a texture on an object (and converts it into a vtkPolyData if necessary).
class VtkCompositeTextureOnSurfaceFilter : public VtkCompositeFilter
{
public:
    VtkCompositeTextureOnSurfaceFilter(vtkAlgorithm* inputAlgorithm);
    virtual ~VtkCompositeTextureOnSurfaceFilter();

    virtual void init() override;

    virtual void SetUserProperty(QString name, QVariant value) override;

private:
};
