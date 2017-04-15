/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-21
 * \brief  Implementation of the VtkCompositeColormapToImageFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeColormapToImageFilter.h"

#include <vtkImageMapToColors.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <QSettings>
#include <QFileDialog>

#include "VtkColorLookupTable.h"
#include "XmlIO/Qt/XmlLutReader.h"

VtkCompositeColormapToImageFilter::VtkCompositeColormapToImageFilter( vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm)
{
    this->init();
}

VtkCompositeColormapToImageFilter::~VtkCompositeColormapToImageFilter() =
    default;

void VtkCompositeColormapToImageFilter::init()
{
    this->_inputDataObjectType = VTK_IMAGE_DATA;
    this->_outputDataObjectType = VTK_IMAGE_DATA;

    vtkSmartPointer<VtkColorLookupTable> colormap = vtkSmartPointer<VtkColorLookupTable>::New();

    QWidget* parent = nullptr;
    QSettings settings;
    QString fileName = QFileDialog::getOpenFileName(parent, "Select color lookup table",
                                                    settings.value("lastOpenedLookupTableFileDirectory").toString(),
                                                    "Lookup table XML files (*.xml);;");
    double range[2];
    dynamic_cast<vtkImageAlgorithm*>(_inputAlgorithm)->GetOutput()->GetPointData()->GetScalars()->GetRange(range);

    DataHolderLib::ColorLookupTable lut;
    if (FileIO::XmlLutReader::readFromFile(fileName, lut))
    {
        colormap->setLookupTable(lut);
        settings.setValue("lastOpenedLookupTableFileDirectory", fileName);
    }
    else
    {
        colormap->SetTableRange(range[0], range[1]);
        colormap->SetHueRange(0.0, 0.666);
    }
    colormap->SetNumberOfTableValues(256);
    colormap->Build();

    colormap->GetTableRange(range);
    QList<QVariant> tableRangeList;
    tableRangeList.push_back(range[0]);
    tableRangeList.push_back(range[1]);
    QList<QVariant> hueRangeList;
    hueRangeList.push_back(0.0);
    hueRangeList.push_back(0.666);
    (*_algorithmUserVectorProperties)["TableRange"] = tableRangeList;
    (*_algorithmUserVectorProperties)["HueRange"] = hueRangeList;

    vtkImageMapToColors* map = vtkImageMapToColors::New();
    map->SetInputConnection(0, _inputAlgorithm->GetOutputPort());
    map->SetLookupTable(colormap);
    map->SetPassAlphaToOutput(1);
    (*_algorithmUserProperties)["PassAlphaToOutput"] = true;
    (*_algorithmUserProperties)["NumberOfColors"] = 256;

    _outputAlgorithm = map;
}

void VtkCompositeColormapToImageFilter::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);

    vtkImageMapToColors* map = static_cast<vtkImageMapToColors*>(_outputAlgorithm);
    if (name.compare("PassAlphaToOutput") == 0)
        map->SetPassAlphaToOutput(value.toBool());
    else if (name.compare("NumberOfColors") == 0)
        static_cast<vtkLookupTable*>(map->GetLookupTable())->SetNumberOfTableValues(value.toInt());
}

void VtkCompositeColormapToImageFilter::SetUserVectorProperty( QString name, QList<QVariant> values )
{
    VtkAlgorithmProperties::SetUserVectorProperty(name, values);

    vtkImageMapToColors* map = static_cast<vtkImageMapToColors*>(_outputAlgorithm);
    if (name.compare("TableRange") == 0)
        static_cast<vtkLookupTable*>(map->GetLookupTable())->SetTableRange(
                values[0].toInt(), values[1].toInt());
    else if (name.compare("HueRange") == 0)
        static_cast<vtkLookupTable*>(map->GetLookupTable())->SetHueRange(values[0].toDouble(),
                                     values[1].toDouble());
}
