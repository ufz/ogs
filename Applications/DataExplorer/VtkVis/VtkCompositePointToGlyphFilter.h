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

#ifndef VTKCOMPOSITEPOINTTOGLYPHFILTER_H
#define VTKCOMPOSITEPOINTTOGLYPHFILTER_H

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief Converts point data to scalar-scaled spheres.
class VtkCompositePointToGlyphFilter : public VtkCompositeFilter
{
public:
    VtkCompositePointToGlyphFilter(vtkAlgorithm* inputAlgorithm);
    virtual ~VtkCompositePointToGlyphFilter();

    virtual void init();

    virtual void SetUserProperty(QString name, QVariant value);

private:
    vtkSphereSource* _glyphSource;
};

#endif // VTKCOMPOSITEPOINTTOGLYPHFILTER_H
