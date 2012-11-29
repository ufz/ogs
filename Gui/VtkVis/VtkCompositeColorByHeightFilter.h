/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkCompositeColorByHeightFilter.h
 *
 * Created on 2010-11-01 by Karsten Rink
 */

#ifndef VTKCOMPOSITECOLORBYHEIGHTFILTER_H
#define VTKCOMPOSITECOLORBYHEIGHTFILTER_H

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief This filter colors the input by the points z-value.
class VtkCompositeColorByHeightFilter : public VtkCompositeFilter
{
public:
	VtkCompositeColorByHeightFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeColorByHeightFilter() {}

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

protected:
};

#endif // VTKCOMPOSITECOLORBYHEIGHTFILTER_H
