/**
 * \file VtkCompositeLineToTubeFilter.h
 * 18/11/2010 KR Initial implementation
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
	int GetInitialRadius() const;
};

#endif // VTKCOMPOSITELINETOTUBEFILTER_H
