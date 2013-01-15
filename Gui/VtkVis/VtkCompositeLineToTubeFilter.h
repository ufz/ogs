/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-18
 * \brief  Definition of the VtkCompositeLineToTubeFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
