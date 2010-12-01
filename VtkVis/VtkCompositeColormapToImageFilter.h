/**
 * \file VtkCompositeColormapToImageFilter.h
 * 21/10/2010 LB Initial implementation
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
