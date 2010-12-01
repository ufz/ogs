/**
 * \file VtkGeoImageSource.h
 * 28/09/2010 LB Initial implementation
 */

#ifndef VTKGEOIMAGESOURCE_H
#define VTKGEOIMAGESOURCE_H

#include <vtkSimpleImageToImageFilter.h>
#include "VtkAlgorithmProperties.h"

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

	std::pair<double, double> getOrigin();

	double getSpacing();

	void setImageFilename(QString filename);
	
	void setImage(QImage& image);
	
	void setOrigin(QPointF origin);
	
	void setSpacing(double spacing);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	/// @brief Constructor.
	VtkGeoImageSource();

	/// @brief Destructor.
	virtual ~VtkGeoImageSource();
	
	/// @brief Filter execution.
	virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);

private:
	VtkGeoImageSource(const VtkGeoImageSource&);	// Not implemented.
	void operator=(const VtkGeoImageSource&);	// Not implemented
	
	vtkQImageToImageSource* _imageSource;
	vtkImageChangeInformation* _imageInfo;
	vtkImageShiftScale* _imageShift;
};

#endif // VTKGEOIMAGESOURCE_H
