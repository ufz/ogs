/**
 * \file VtkVisHelper.h
 * 22/09/2010 LB Initial implementation
 */

#ifndef VTKVISHELPER_H
#define VTKVISHELPER_H

class QImage;
class vtkTexture;
class vtkImageData;

/**
 * \brief Some data conversion helper functions.
 */
class VtkVisHelper
{
public:
	/// @brief Converts a QImage to vtkImageData.
	static vtkImageData* QImageToVtkImageData(QImage &img);

	/// @brief Converts a QImage-object into a vtkTexture-object
	static vtkTexture* QImageToVtkTexture(QImage &img);
};

#endif // VTKVISHELPER_H
