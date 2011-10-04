/**
 * \file VtkVisImageItem.cpp
 * 2011/09/29 KR Initial implementation
 */

// ** INCLUDES **
#include "VtkVisImageItem.h"
#include "VtkAlgorithmProperties.h"

#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include "QVtkDataSetMapper.h"
#include <vtkImageAlgorithm.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

// export test
#include <vtkXMLImageDataWriter.h>
#include <vtkImageActor.h>


VtkVisImageItem::VtkVisImageItem(
	vtkAlgorithm* algorithm, TreeItem* parentItem,
	const QList<QVariant> data /*= QList<QVariant>()*/)
: VtkVisPipelineItem(algorithm, parentItem, data)
{
}

VtkVisImageItem::VtkVisImageItem(
	VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
	const QList<QVariant> data /*= QList<QVariant>()*/)
: VtkVisPipelineItem(compositeFilter, parentItem, data)
{
}

VtkVisImageItem::~VtkVisImageItem()
{
}

void VtkVisImageItem::Initialize(vtkRenderer* renderer)
{
	_renderer = renderer;
	_mapper = QVtkDataSetMapper::New();
	_mapper->InterpolateScalarsBeforeMappingOff();

	// Use a special vtkImageActor instead of vtkActor
	vtkImageAlgorithm* imageAlg = static_cast<vtkImageAlgorithm*>(_algorithm);
	vtkImageActor* imageActor = vtkImageActor::New();
	imageActor->SetInput(imageAlg->GetOutput());
	_actor = imageActor;
	_renderer->AddActor(_actor);

	// Set pre-set properties
	VtkAlgorithmProperties* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(_algorithm);
	if (vtkProps)
		setVtkProperties(vtkProps);
/*
	// Copy properties from parent
	else
	{
*/
		VtkVisPipelineItem* parentItem = dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
		while (parentItem)
		{
			VtkAlgorithmProperties* parentProps = dynamic_cast<VtkAlgorithmProperties*>(parentItem->algorithm());
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
//	}

	// Set active scalar to the desired one from VtkAlgorithmProperties
	// or to match those of the parent.
	if (vtkProps)
	{
		if (vtkProps->GetActiveAttribute().length() > 0)
		{
			this->SetActiveAttribute(vtkProps->GetActiveAttribute());
		}
		else
		{
			VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
			if (visParentItem)
				this->SetActiveAttribute(visParentItem->GetActiveAttribute());
			if (vtkProps->GetTexture() != NULL)
				this->SetActiveAttribute("Solid Color");
		}
	}
}

void VtkVisImageItem::setVtkProperties(VtkAlgorithmProperties* vtkProps)
{
	// todo
}

const int VtkVisImageItem::callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const
{
	vtkImageAlgorithm* algID = dynamic_cast<vtkImageAlgorithm*>(algorithm);
	if (algID)
	{
		vtkSmartPointer<vtkXMLImageDataWriter> iWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
		iWriter->SetInput(algID->GetOutputDataObject(0));
		std::string filenameWithExt = filename;
		filenameWithExt.append(".vti");
		iWriter->SetFileName(filenameWithExt.c_str());
		return iWriter->Write();
	}
	std::cout << "VtkVisPipelineItem::writeToFile() - Unknown data type..." << std::endl;
	return 0;
}
