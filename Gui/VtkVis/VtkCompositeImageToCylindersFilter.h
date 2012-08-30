/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VtkCompositeImageToCylindersFilter.h
 *
 * Created on 2010-10-19 by Lars Bilke
 */

#ifndef VTKCOMPOSITEIMAGETOCYLINDERSFILTER_H
#define VTKCOMPOSITEIMAGETOCYLINDERSFILTER_H

#include "VtkCompositeFilter.h"

class VtkImageDataToLinePolyDataFilter;

/// @brief Creates cylinders that stand on top of the image with the length
/// of the corresponding first sub-pixel value (the red value). Useful to
/// visualize precipitation maps as a 3d bar chart.
class VtkCompositeImageToCylindersFilter : public VtkCompositeFilter
{
public:
	VtkCompositeImageToCylindersFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeImageToCylindersFilter();

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

	void SetUserVectorProperty( QString name, QList<QVariant> values );

private:
	VtkImageDataToLinePolyDataFilter* _lineFilter;
};

#endif // VTKCOMPOSITEIMAGETOCYLINDERSFILTER_H
