/**
 * \file VtkImageDataToLinePolyDataFilter.h
 * 06/10/2010 LB Initial implementation
 */

#ifndef VTKIMAGEDATATOLINEPOLYDATAFILTER_H
#define VTKIMAGEDATATOLINEPOLYDATAFILTER_H

#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

class vtkInformation;
class vtkInformationVector;

/// @brief Creates lines that stand on top of the image with the length
/// of the corresponding first sub-pixel value (the grey or red value).
/// The maximum height is 0.1 * longest image dimension.
/// Used by VtkCompositeImageDataToCylindersFilter.
class VtkImageDataToLinePolyDataFilter : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
	/// @brief Create new objects with New() because of VTKs reference counting.
	static VtkImageDataToLinePolyDataFilter* New();

	vtkTypeMacro(VtkImageDataToLinePolyDataFilter, vtkPolyDataAlgorithm);

	/// @brief Prints information about itself.
	void PrintSelf(ostream& os, vtkIndent indent);

	/// @brief Sets the scaling of the length of the lines.
	ogsUserPropertyMacro(LengthScaleFactor,double);

	/// @brief Sets a user property.
	virtual void SetUserProperty(QString name, QVariant value)
	{
		if (name.compare("LengthScaleFactor") == 0)
			SetLengthScaleFactor(value.toDouble());
	}

	/// @brief Returns the space between two pixels.
	vtkGetMacro(ImageSpacing,double);

protected:
	/// @brief Constructor.
	VtkImageDataToLinePolyDataFilter();

	/// @brief Destructor.
	virtual ~VtkImageDataToLinePolyDataFilter();

	/// @brief Sets input port to vtkImageData.
	virtual int FillInputPortInformation(int port, vtkInformation* info);

	/// @brief Converts the image data to lines
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
	                        vtkInformationVector* outputVector);

	/// @brief The spacing of the image
	double ImageSpacing;

private:
	VtkImageDataToLinePolyDataFilter(const VtkImageDataToLinePolyDataFilter&); // Not implemented.
	void operator=(const VtkImageDataToLinePolyDataFilter&); // Not implemented
};

#endif // VTKIMAGEDATATOLINEPOLYDATAFILTER_H
