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


class QImage;

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

	/// Sets the surfaces vector
	void SetImage(vtkTexture* texture);

	/// Sets the cellsize (i.e. the actual dimension of a pixel)
	void SetCellSize(double c) { _cellsize = c; };

	/// Sets the raster/image to be used as a texture map
	void SetRaster(QImage &img);

	/// Sets the geo-referenced origin of the image (i.e. the lower left corner)
	virtual void SetOrigin(double x, double y, double z = 0.0) { _origin.first = x; _origin.second = y; (void)z; };
	virtual void SetOrigin(double* pos) { _origin.first = pos[0]; _origin.second = pos[1]; };

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkBGImageSource();
	~VtkBGImageSource();


private:
	
	std::pair<double, double> _origin;
	double _cellsize;


};

#endif // VTKBGIMAGESOURCE_H
