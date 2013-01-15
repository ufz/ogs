/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-21
 * \brief  Definition of the VtkCompositeColormapToImageFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
