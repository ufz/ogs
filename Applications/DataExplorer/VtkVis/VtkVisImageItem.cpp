/**
 * \file
 * \author Karsten Rink
 * \date   2011-09-29
 * \brief  Implementation of the VtkVisImageItem class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisImageItem.h"

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"

#include "VtkAlgorithmProperties.h"
#include "VtkGeoImageSource.h"

#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkImageAlgorithm.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkImageShiftScale.h>
#include <vtkDataArray.h>

// export
#include <vtkImageActor.h>
#include <vtkXMLImageDataWriter.h>

VtkVisImageItem::VtkVisImageItem(
    vtkAlgorithm* algorithm, TreeItem* parentItem,
    const QList<QVariant> data /*= QList<QVariant>()*/)
    : VtkVisPipelineItem(algorithm, parentItem, data), _transformFilter(nullptr)
{
}

VtkVisImageItem::VtkVisImageItem(
    VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
    const QList<QVariant> data /*= QList<QVariant>()*/)
    : VtkVisPipelineItem(compositeFilter, parentItem, data),
      _transformFilter(nullptr)
{
}

VtkVisImageItem::~VtkVisImageItem()
{
    _transformFilter->Delete();
}

void VtkVisImageItem::Initialize(vtkRenderer* renderer)
{
    auto* img = dynamic_cast<vtkImageAlgorithm*>(_algorithm);
    img->Update();
    // VtkGeoImageSource* img = dynamic_cast<VtkGeoImageSource*>(_algorithm);

    double origin[3];
    double spacing[3];
    double range[2];
    img->GetOutput()->GetOrigin(origin);
    img->GetOutput()->GetSpacing(spacing);
    // img->getRange(range);
    img->GetOutput()->GetPointData()->GetScalars()->GetRange(range);
    vtkImageShiftScale* scale = vtkImageShiftScale::New();
    scale->SetOutputScalarTypeToUnsignedChar();
    scale->SetInputConnection(img->GetOutputPort());
    scale->SetShift(-range[0]);
    scale->SetScale(255.0 / (range[1] - range[0]));

    _transformFilter = vtkImageChangeInformation::New();
    _transformFilter->SetInputConnection(scale->GetOutputPort());
    // double origin[3];
    // img->getOrigin(origin);
    // double spacing = img->getSpacing();
    //_transformFilter->SetOutputOrigin(origin2);
    //_transformFilter->SetOutputSpacing(spacing2);
    _transformFilter->Update();

    _renderer = renderer;

    // Use a special vtkImageActor instead of vtkActor
    vtkImageActor* imageActor = vtkImageActor::New();
    imageActor->SetInputData(_transformFilter->GetOutput());
    _actor = imageActor;
    _renderer->AddActor(_actor);

    // Set pre-set properties
    auto* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(_algorithm);
    if (vtkProps)
        setVtkProperties(vtkProps);

    auto* parentItem = dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
    while (parentItem)
    {
        auto* parentProps =
            dynamic_cast<VtkAlgorithmProperties*>(parentItem->algorithm());
        if (parentProps)
        {
            auto* newProps = new VtkAlgorithmProperties();
            newProps->SetScalarVisibility(parentProps->GetScalarVisibility());
            newProps->SetTexture(parentProps->GetTexture());
            setVtkProperties(newProps);
            vtkProps = newProps;
            parentItem = nullptr;
        }
        else
            parentItem =
                dynamic_cast<VtkVisPipelineItem*>(parentItem->parentItem());
    }

    // Set active scalar to the desired one from VtkAlgorithmProperties
    // or to match those of the parent.
    if (vtkProps)
    {
        if (vtkProps->GetActiveAttribute().length() > 0)
            this->SetActiveAttribute(vtkProps->GetActiveAttribute());
        else
        {
            auto* visParentItem =
                dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
            if (visParentItem)
                this->SetActiveAttribute(visParentItem->GetActiveAttribute());
            if (vtkProps->GetTexture() != nullptr)
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
    std::string file_name_cpy(filename);
    auto* algID = dynamic_cast<vtkImageAlgorithm*>(algorithm);
    if (algID)
    {
        vtkSmartPointer<vtkXMLImageDataWriter> iWriter =
                vtkSmartPointer<vtkXMLImageDataWriter>::New();
        iWriter->SetInputData(algID->GetOutputDataObject(0));
        if (BaseLib::getFileExtension(filename) != "vti")
            file_name_cpy.append(".vti");
        iWriter->SetFileName(file_name_cpy.c_str());
        return iWriter->Write();
    }
    ERR("VtkVisPipelineItem::writeToFile() - Unknown data type.");
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
