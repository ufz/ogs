/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VtkCompositeLineToTubeFilter.h
 *
 * Created on 2010-11-18 by Karsten Rink
 */

#ifndef VTKCOMPOSITELINETOTUBEFILTER_H
#define VTKCOMPOSITELINETOTUBEFILTER_H

#include "VtkCompositeFilter.h"

/// @brief Converts lines to tube-objects.
class VtkCompositeLineToTubeFilter : public VtkCompositeFilter
{
public:
	VtkCompositeLineToTubeFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeLineToTubeFilter();

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

private:
	float GetInitialRadius() const;
};

#endif // VTKCOMPOSITELINETOTUBEFILTER_H
