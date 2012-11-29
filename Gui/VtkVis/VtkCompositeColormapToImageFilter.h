/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkCompositeColormapToImageFilter.h
 *
 * Created on 2010-10-21 by Lars Bilke
 */

#ifndef VTKCOMPOSITECOLORMAPTOIMAGEFILTER_H
#define VTKCOMPOSITECOLORMAPTOIMAGEFILTER_H

#include "VtkCompositeFilter.h"

/// @brief Applies a user adjustable color map to an image.
class VtkCompositeColormapToImageFilter : public VtkCompositeFilter
{
public:
	VtkCompositeColormapToImageFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeColormapToImageFilter();

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

	virtual void SetUserVectorProperty(QString name, QList<QVariant> values);

private:
};

#endif // VTKCOMPOSITECOLORMAPTOIMAGEFILTER_H
