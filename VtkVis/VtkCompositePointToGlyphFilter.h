/**
 * \file VtkCompositePointToGlyphFilter.h
 * 21/10/2010 LB Initial implementation
 */

#ifndef VTKCOMPOSITEPOINTTOGLYPHFILTER_H
#define VTKCOMPOSITEPOINTTOGLYPHFILTER_H

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief Converts point data to scalar-scaled spheres.
class VtkCompositePointToGlyphFilter : public VtkCompositeFilter
{
public:
	VtkCompositePointToGlyphFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositePointToGlyphFilter();

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

private:
	vtkSphereSource* _glyphSource;
	
};

#endif // VTKCOMPOSITEPOINTTOGLYPHFILTER_H
