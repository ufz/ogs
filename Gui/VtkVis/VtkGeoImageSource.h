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

/**
 * \brief The VtkVisPipeline source object of a geo-referenced image (file).
 */
class VtkGeoImageSource : public vtkSimpleImageToImageFilter, public VtkAlgorithmProperties
{
public:
	/// @brief Create new objects with New() because of VTKs reference counting.
	static VtkGeoImageSource* New();

	vtkTypeMacro(VtkGeoImageSource, vtkSimpleImageToImageFilter);

	/// @brief Prints information about itself.
	void PrintSelf(ostream& os, vtkIndent indent);

	/// @brief Returns the ImageData object.
	vtkImageData* getImageData();

	/// @brief Reads an image from file.
	void readImage(const QString &filename);

	/// @brief Imports an existing image object.
	void setImage(vtkImageAlgorithm* img, const QString &name, double x0, double y0, double spacing);

	/// @brief Returns the origin in world coordinates.
	void getOrigin(double origin[3]) const;

	/// @brief Returns the scalar data range.
	void getRange(double range[2]);

	/// @brief Returns the spacing betweeen two pixels.
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

	double _x0, _y0, _z0, _spacing;
};

#endif // VTKGEOIMAGESOURCE_H
