/**
 * \file VtkCompositeTextureOnSurfaceFilter.h
 * 18/11/2010 KR Initial implementation
 */

#ifndef VTKCOMPOSITETEXTUREONSURFACEFILTER_H
#define VTKCOMPOSITETEXTUREONSURFACEFILTER_H

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief Puts a texture on an object (and converts it into a vtkPolyData if necessary).
class VtkCompositeTextureOnSurfaceFilter : public VtkCompositeFilter
{
public:
	VtkCompositeTextureOnSurfaceFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeTextureOnSurfaceFilter();

	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

private:
};

#endif // VTKCOMPOSITETEXTUREONSURFACEFILTER_H
