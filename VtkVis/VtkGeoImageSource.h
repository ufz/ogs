/**
 * \file VtkGeoImageSource.h
 * 28/09/2010 LB Initial implementation
 */

#ifndef VTKGEOIMAGESOURCE_H
#define VTKGEOIMAGESOURCE_H

#include "VtkAlgorithmProperties.h"
#include <vtkSimpleImageToImageFilter.h>

class QString;
class QPointF;
class QImage;
class vtkQImageToImageSource;
class vtkImageShiftScale;
class vtkImageData;

class VtkGeoImageSource : public vtkSimpleImageToImageFilter, public VtkAlgorithmProperties
{
public:
	/// @brief Create new objects with New() because of VTKs reference counting.
	static VtkGeoImageSource* New();

	vtkTypeMacro(VtkGeoImageSource, vtkSimpleImageToImageFilter);

	/// @brief Prints information about itself.
	void PrintSelf(ostream& os, vtkIndent indent);

	vtkImageData* getImageData();

	void setImageFilename(QString filename);

	void getOrigin(double& x, double& y, double& z) const;

	double getSpacing() const;

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	/// @brief Constructor.
	VtkGeoImageSource();

	/// @brief Destructor.
	virtual ~VtkGeoImageSource();

	/// @brief Filter execution.
	virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);

private:
	VtkGeoImageSource(const VtkGeoImageSource&); // Not implemented.
	void operator=(const VtkGeoImageSource&); // Not implemented

	vtkImageAlgorithm* _imageSource;
	vtkImageShiftScale* _imageShift;

	double _x0, _y0, _z0, _spacing;
};

#endif // VTKGEOIMAGESOURCE_H
