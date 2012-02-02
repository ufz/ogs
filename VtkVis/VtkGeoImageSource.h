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
class vtkImageChangeInformation;
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

	vtkImageChangeInformation* getImageChangeTransformation() const { return _imageInfo; };

	void setImageFilename(QString filename);

	//void setImage(QImage& image);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	/// @brief Constructor.
	VtkGeoImageSource();

	/// @brief Destructor.
	virtual ~VtkGeoImageSource();

	/// @brief Filter execution.
	virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);

	void setOrigin(double x, double y, double z);

	void setSpacing(double spacing);


private:
	VtkGeoImageSource(const VtkGeoImageSource&); // Not implemented.
	void operator=(const VtkGeoImageSource&); // Not implemented

	vtkImageAlgorithm* _imageSource;
	vtkImageChangeInformation* _imageInfo;
	vtkImageShiftScale* _imageShift;
};

#endif // VTKGEOIMAGESOURCE_H
