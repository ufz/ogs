/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkCompositePointToGlyphFilter.h
 *
 * Created on 2010-10-21 by Lars Bilke
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
	float GetInitialRadius() const;

	vtkSphereSource* _glyphSource;
};

#endif // VTKCOMPOSITEPOINTTOGLYPHFILTER_H
