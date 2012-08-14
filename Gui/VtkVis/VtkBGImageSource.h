/**
 * \file VtkBGImageSource.h
 * 30/04/2010 KR Initial implementation
 *
 */

#ifndef VTKBGIMAGESOURCE_H
#define VTKBGIMAGESOURCE_H

// ** INCLUDES **
#include <vtkTextureMapToPlane.h>

#include "VtkAlgorithmProperties.h"

class vtkImageAlgorithm;

/**
 * \brief Uses an image source to create a plane in the 3D with the given
 * image texture mapped on it.
 */
class VtkBGImageSource : public vtkTextureMapToPlane, public VtkAlgorithmProperties
{
public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkBGImageSource* New();

	vtkTypeRevisionMacro(VtkBGImageSource, vtkTextureMapToPlane);

	/// Sets the raster/image to be used as a texture map
	void SetRaster(vtkImageAlgorithm *img, double x0, double y0, double scalingFactor);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkBGImageSource();
	~VtkBGImageSource();

private:

	std::pair<double, double> _origin;
	double _cellsize;
};

#endif // VTKBGIMAGESOURCE_H
