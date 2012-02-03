/**
 * \file VtkVisImageItem.cpp
 * 2011/09/29 KR Initial implementation
 */

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include "VtkVisImageItem.h"
#include "VtkGeoImageSource.h"

#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkImageAlgorithm.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

// export
#include <vtkImageActor.h>
#include <vtkXMLImageDataWriter.h>

VtkVisImageItem::VtkVisImageItem(
        vtkAlgorithm* algorithm, TreeItem* parentItem,
        const QList<QVariant> data /*= QList<QVariant>()*/)
	: VtkVisPipelineItem(algorithm, parentItem, data), _transformFilter(NULL)
{
}

VtkVisImageItem::VtkVisImageItem(
        VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
        const QList<QVariant> data /*= QList<QVariant>()*/)
	: VtkVisPipelineItem(compositeFilter, parentItem, data), _transformFilter(NULL)
{
}

VtkVisImageItem::~VtkVisImageItem()
{
	_transformFilter->Delete();
}

void VtkVisImageItem::Initialize(vtkRenderer* renderer)
{
	VtkGeoImageSource* img = dynamic_cast<VtkGeoImageSource*>(_algorithm);

	_transformFilter = vtkImageChangeInformation::New();
	_transformFilter->SetInputConnection(img->GetOutputPort());
	double x0, y0, z0;
	img->getOrigin(x0, y0, z0);
	double spacing = img->getSpacing();
	_transformFilter->SetOutputOrigin(x0, y0, z0);
	_transformFilter->SetOutputSpacing(spacing, spacing, spacing);
	_transformFilter->Update();

	_renderer = renderer;

	// Use a special vtkImageActor instead of vtkActor
	vtkImageActor* imageActor = vtkImageActor::New();
	imageActor->SetInput(_transformFilter->GetOutput());
	_actor = imageActor;
	_renderer->AddActor(_actor);

	// Set pre-set properties
	VtkAlgorithmProperties* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(_algorithm);
	if (vtkProps)
		setVtkProperties(vtkProps);

	VtkVisPipelineItem* parentItem = dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
	while (parentItem)
	{
		VtkAlgorithmProperties* parentProps =
		        dynamic_cast<VtkAlgorithmProperties*>(parentItem->algorithm());
		if (parentProps)
		{
			VtkAlgorithmProperties* newProps = new VtkAlgorithmProperties();
			newProps->SetScalarVisibility(parentProps->GetScalarVisibility());
			newProps->SetTexture(parentProps->GetTexture());
			setVtkProperties(newProps);
			vtkProps = newProps;
			parentItem = NULL;
		}
		else
			parentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem->parentItem());
	}

	// Set active scalar to the desired one from VtkAlgorithmProperties
	// or to match those of the parent.
	if (vtkProps)
	{
		if (vtkProps->GetActiveAttribute().length() > 0)
			this->SetActiveAttribute(vtkProps->GetActiveAttribute());
		else
		{
			VtkVisPipelineItem* visParentItem =
			        dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
			if (visParentItem)
				this->SetActiveAttribute(visParentItem->GetActiveAttribute());
			if (vtkProps->GetTexture() != NULL)
				this->SetActiveAttribute("Solid Color");
		}
	}
}

void VtkVisImageItem::setVtkProperties(VtkAlgorithmProperties* vtkProps)
{
	// todo implementation
	(void)vtkProps;
}

int VtkVisImageItem::callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const
{
	vtkImageAlgorithm* algID = dynamic_cast<vtkImageAlgorithm*>(algorithm);
	if (algID)
	{
		vtkSmartPointer<vtkXMLImageDataWriter> iWriter =
		        vtkSmartPointer<vtkXMLImageDataWriter>::New();
		iWriter->SetInput(algID->GetOutputDataObject(0));
		std::string filenameWithExt = filename;
		filenameWithExt.append(".vti");
		iWriter->SetFileName(filenameWithExt.c_str());
		return iWriter->Write();
	}
	std::cout << "VtkVisPipelineItem::writeToFile() - Unknown data type..." << std::endl;
	return 0;
}

void VtkVisImageItem::setTranslation(double x, double y, double z) const
{
	_transformFilter->SetOriginTranslation(x,y,z);
}

vtkAlgorithm* VtkVisImageItem::transformFilter() const
{
	return _transformFilter;
}
