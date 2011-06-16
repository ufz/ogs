/**
 * \file VtkTextureOnSurfaceFilter.cpp
 * 3/2/2010 LB Initial implementation
 * 23/04/2010 KR Surface visualisation
 *
 * Implementation of VtkSurfacesSource
 */

// ** INCLUDES **
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>

#include "MathTools.h"
#include "VtkTextureOnSurfaceFilter.h"
#include "VtkVisHelper.h"

#include <QImage>


vtkStandardNewMacro(VtkTextureOnSurfaceFilter);
vtkCxxRevisionMacro(VtkTextureOnSurfaceFilter, "$Revision$");

VtkTextureOnSurfaceFilter::VtkTextureOnSurfaceFilter() : _origin(0,0), _scalingFactor(0)
{
}

VtkTextureOnSurfaceFilter::~VtkTextureOnSurfaceFilter()
{
}

void VtkTextureOnSurfaceFilter::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);
}

int VtkTextureOnSurfaceFilter::RequestData( vtkInformation* request,
							             vtkInformationVector** inputVector,
								         vtkInformationVector* outputVector )
{
	(void)request;

	if (this->GetTexture() == NULL)
	{
		std::cout << "Error in VtkTextureOnSurfaceFilter::RequestData() - No texture specified ..." << std::endl;
		return 0;
	}

	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	int dims[3];
	this->GetTexture()->GetInput()->GetDimensions(dims);
	size_t imgWidth=dims[0]; size_t imgHeight=dims[1];

	std::pair<int, int> min((int)_origin.first, (int)_origin.second);
	std::pair<int, int> max((int)(_origin.first+(imgWidth * _scalingFactor)), (int)(_origin.second+(imgHeight * _scalingFactor)));

	//calculate texture coordinates
	vtkPoints* points = input->GetPoints();
	vtkSmartPointer<vtkFloatArray> textureCoordinates = vtkSmartPointer<vtkFloatArray>::New();
		textureCoordinates->SetNumberOfComponents(3);
		textureCoordinates->SetName("TextureCoordinates");
		size_t nPoints = points->GetNumberOfPoints();
		for (size_t i=0; i<nPoints; i++)
		{
			double coords[3];
			points->GetPoint(i, coords);
			float newcoords[3] = {MathLib::normalize(min.first, max.first, coords[0]), MathLib::normalize(min.second, max.second, coords[1]), 0 /*coords[2]*/ };
			textureCoordinates->InsertNextTuple(newcoords);
		}

	// put it all together
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
		output->CopyStructure(input);
		output->GetPointData()->PassData(input->GetPointData());
		output->GetCellData()->PassData(input->GetCellData());
		output->GetPointData()->SetTCoords(textureCoordinates);


	return 1;
}

void VtkTextureOnSurfaceFilter::SetRaster(QImage &img, std::pair<float, float> origin, double scalingFactor)
{
	_origin = origin;
	_scalingFactor = scalingFactor;
	QImage raster = img.transformed(QTransform(1, 0, 0, -1, 0, 0), Qt::FastTransformation);
	this->SetTexture(VtkVisHelper::QImageToVtkTexture(raster));
}

void VtkTextureOnSurfaceFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);
}
