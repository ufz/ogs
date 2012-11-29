/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkCompositeThresholdFilter.h
 *
 * Created on 2010-10-15 by Lars Bilke
 */

#ifndef VTKCOMPOSITETHRESHOLDFILTER_H
#define VTKCOMPOSITETHRESHOLDFILTER_H

#include "VtkCompositeFilter.h"

/// @brief Visualises only parts of meshes that are above/below/within given thresholds.
/// In init() the threshold is first set to double min / max values. Set the
/// threshold later on via SetUserVectorProperty() to the actual data range.
class VtkCompositeThresholdFilter : public VtkCompositeFilter
{
public:
	VtkCompositeThresholdFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeThresholdFilter();

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

	void SetUserVectorProperty( QString name, QList<QVariant> values );

private:
};

#endif // VTKCOMPOSITETHRESHOLDFILTER_H
