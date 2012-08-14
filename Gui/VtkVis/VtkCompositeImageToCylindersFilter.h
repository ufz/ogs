/**
 * \file VtkCompositeImageToCylindersFilter.h
 * 19/10/2010 LB Initial implementation
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
