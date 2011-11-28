/**
 * \file VtkGeoImageSource.cpp
 * 28/09/2010 LB Initial implementation
 *
 * Implementation of VtkGeoImageSource class
 */

// ** INCLUDES **
#include "VtkGeoImageSource.h"

#include "OGSRaster.h"

#include <vtkFloatArray.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkQImageToImageSource.h>

#include <QImage>
#include <QPointF>
#include <QString>

vtkStandardNewMacro(VtkGeoImageSource);

// This function is copied from the vtkSimpleImageFilterExample.cpp
//
// The switch statement in Execute will call this method with
// the appropriate input type (IT). Note that this example assumes
// that the output data type is the same as the input data type.
// This is not always the case.
template <class IT>
void vtkSimpleImageFilterExampleExecute(vtkImageData* input,
                                        vtkImageData* output,
                                        IT* inPtr, IT* outPtr)
{
	int dims[3];
	input->GetDimensions(dims);
	if (input->GetScalarType() != output->GetScalarType())
	{
		vtkGenericWarningMacro(<< "Execute: input ScalarType, " << input->GetScalarType()
		                       << ", must match out ScalarType " << output->GetScalarType());
		return;
	}
	// HACK LB Multiply by number of scalar components due to RGBA values ?????
	int size = dims[0] * dims[1] * dims[2] * input->GetNumberOfScalarComponents();

	for(int i = 0; i < size; i++)
		outPtr[i] = inPtr[i];
}

VtkGeoImageSource::VtkGeoImageSource()
{
	_imageSource = vtkQImageToImageSource::New();
	_imageInfo = vtkImageChangeInformation::New();
	_imageShift = vtkImageShiftScale::New();

	_imageInfo->SetInputConnection(_imageSource->GetOutputPort());
	_imageShift->SetInputConnection(_imageInfo->GetOutputPort());
	_imageShift->SetOutputScalarTypeToUnsignedChar();
	this->SetInputConnection(_imageShift->GetOutputPort());
}

VtkGeoImageSource::~VtkGeoImageSource()
{
	_imageSource->Delete();
	_imageInfo->Delete();
	_imageShift->Delete();
}

void VtkGeoImageSource::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

void VtkGeoImageSource::setImageFilename(QString filename)
{
	QImage* raster = new QImage;
	QPointF origin;
	double spacing;
	OGSRaster::loadImage(filename, *raster, origin, spacing);
	this->setImage(*raster);
	delete raster;
	// correct raster position by half a pixel for correct visualisation 
	this->setOrigin(origin.x()+(spacing/2.0), origin.y()+(spacing/2.0), -10.0);
	this->setSpacing(spacing);
	this->SetName(filename);
}

vtkImageData* VtkGeoImageSource::getImageData()
{
	return this->_imageSource->GetImageDataInput(0);
}

std::pair<double, double> VtkGeoImageSource::getOrigin()
{
	double* origin = this->_imageInfo->GetOutputOrigin();
	std::pair<double, double> p(origin[0], origin[1]);
	return p;
}

double VtkGeoImageSource::getSpacing()
{
	double* spacing = this->_imageInfo->GetOutputSpacing();
	return spacing[0];
}

void VtkGeoImageSource::setImage(QImage& image)
{
	_imageSource->SetQImage(&image);
	_imageSource->Update(); // crashes otherwise
}

void VtkGeoImageSource::setOrigin(double x, double y, double z)
{
	_imageInfo->SetOutputOrigin(x, y, z);
}

void VtkGeoImageSource::setSpacing(double spacing)
{
	_imageInfo->SetOutputSpacing(spacing, spacing, spacing);
}

void VtkGeoImageSource::SimpleExecute(vtkImageData* input, vtkImageData* output)
{
	vtkDebugMacro(<< "Executing VtkGeoImageSource")
	void* inPtr = input->GetScalarPointer();
	void* outPtr = output->GetScalarPointer();
	switch(output->GetScalarType())
	{
		// This is simply a #define for a big case list.
		// It handles all data types that VTK supports.
		vtkTemplateMacro(vtkSimpleImageFilterExampleExecute(
		                         input, output, (VTK_TT*)(inPtr), (VTK_TT*)(outPtr)));
	default:
		vtkGenericWarningMacro("Execute: Unknown input ScalarType");
		return;
	}

/*
    // Create normals
    vtkFloatArray* normals = vtkFloatArray::New();
    int numPoints = input->GetNumberOfPoints();
    normals->SetNumberOfComponents(3);
    normals->Allocate(3*numPoints);
    float normal[3] = {0.f, 0.f, 1.f};

    // vector data
    std::cout << input->GetScalarTypeAsString() << std::endl;
    vtkFloatArray* vectors = vtkFloatArray::New();
    vectors->SetNumberOfComponents(3);
    vectors->Allocate(3*numPoints);
    int numScalarComponents = input->GetNumberOfScalarComponents();

    for (int i = 0; i < numPoints; i++)
    {
        normals->InsertTuple(i, normal);
        float vector[3];
        vector[0] = 1;
        vector[1] = 1;
        vector[2] = ((unsigned char*)inPtr)[i * numScalarComponents] * 0.1;
        //std::cout << vector[2] << " ";
        vectors->InsertTuple(i, vector);
    }

    normals->SetName("Normals");
    output->GetPointData()->SetNormals(normals);
    normals->Delete();

    vectors->SetName("Vectors");
    output->GetPointData()->SetVectors(vectors);
    vectors->Delete();
 */
}

void VtkGeoImageSource::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}
