/**
 * \file VtkOGSPolyDataAlgorithm.h
 * 28/05/2010 KR initial implementation
 *
 */


#ifndef VTKOGSPOLYDATAALGORITHM_H
#define VTKOGSPOLYDATAALGORITHM_H

// ** INCLUDES **
#include <vtkPolyDataAlgorithm.h>
#include <vtkImageData.h>
#include "VtkAlgorithmProperties.h"

class QImage;

/**
 * \brief Filter class for assigning a texture to a surface.
 *
 * Use SetRaster() to define the texture that should be mapped on the object.
 * The input of this class is a vtkPolyData object. It is important to call SetTexture() before
 * calling SetInputConnection(). Texture coordinates will then be calculated automatically and
 * the texture will also be saved in the VtkAlgorithmProperties object from which this class is
 * derived (i.e. the texture can be returned by VtkAlgorithmProperties::GetTexture()).
 *
 * For convenience this class also has a converter function ConvertImageToTexture() which uses 
 * a QImage as input.
 */
class VtkTextureOnSurfaceFilter : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkTextureOnSurfaceFilter* New();

	vtkTypeRevisionMacro(VtkTextureOnSurfaceFilter,vtkPolyDataAlgorithm);

	/// Prints the object data to an output stream.
	void PrintSelf(ostream& os, vtkIndent indent);	

	/// Sets the raster/image to be used as a texture map
	void SetRaster(QImage &img, std::pair<float, float> origin, double scalingFactor);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkTextureOnSurfaceFilter();
	~VtkTextureOnSurfaceFilter();

	/// Copies the input data and calculates texture coordinates (this requires that SetRaster() has
	/// been called before this method is executed.
	int RequestData(vtkInformation* request, 
		            vtkInformationVector** inputVector, 
					vtkInformationVector* outputVector);

private:
	std::pair<float, float> _origin;
	double _scalingFactor;
};

#endif // VTKOGSPOLYDATAALGORITHM_H
