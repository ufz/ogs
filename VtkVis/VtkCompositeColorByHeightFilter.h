/**
 * \file VtkCompositeColorByHeightFilter.h
 * 01/11/2010 KR Initial implementation
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
	virtual ~VtkCompositeColorByHeightFilter() {};

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	
};

#endif // VTKCOMPOSITECOLORBYHEIGHTFILTER_H
