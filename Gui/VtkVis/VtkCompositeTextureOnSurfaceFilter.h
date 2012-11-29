/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkCompositeTextureOnSurfaceFilter.h
 *
 * Created on 2010-11-18 by Karsten Rink
 */

#ifndef VTKCOMPOSITETEXTUREONSURFACEFILTER_H
#define VTKCOMPOSITETEXTUREONSURFACEFILTER_H

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief Puts a texture on an object (and converts it into a vtkPolyData if necessary).
class VtkCompositeTextureOnSurfaceFilter : public VtkCompositeFilter
{
public:
	VtkCompositeTextureOnSurfaceFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeTextureOnSurfaceFilter();

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

private:
};

#endif // VTKCOMPOSITETEXTUREONSURFACEFILTER_H
