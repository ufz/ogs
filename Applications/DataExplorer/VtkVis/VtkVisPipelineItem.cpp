/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-17
 * \brief  Implementation of the VtkVisPipelineItem class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisPipelineItem.h"

#include <vtkActor.h>
#include <vtkAlgorithm.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkImageActor.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTextureMapToPlane.h>

#include <QMessageBox>

#include "BaseLib/FileTools.h"
#include "QVtkDataSetMapper.h"
#include "VtkAlgorithmProperties.h"
#include "VtkCompositeFilter.h"

VtkVisPipelineItem::VtkVisPipelineItem(
    vtkAlgorithm* algorithm, TreeItem* parentItem,
    const QList<QVariant> data /*= QList<QVariant>()*/)
    : TreeItem(data, parentItem),
      _actor(nullptr),
      _algorithm(algorithm),
      _renderer(nullptr),
      _compositeFilter(nullptr),
      _vtkProps(nullptr)
{
    auto* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
    if (parentItem->parentItem())
    {
        _algorithm->SetInputConnection(
            visParentItem->algorithm()->GetOutputPort());
    }
}

VtkVisPipelineItem::VtkVisPipelineItem(
    VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
    const QList<QVariant> data /*= QList<QVariant>()*/)
    : TreeItem(data, parentItem),
      _actor(nullptr),
      _renderer(nullptr),
      _compositeFilter(compositeFilter),
      _vtkProps(nullptr)
{
    _algorithm = _compositeFilter->GetOutputAlgorithm();
}

VtkVisPipelineItem::~VtkVisPipelineItem()
{
    _renderer->RemoveActor(_actor);
    _actor->Delete();
    delete _compositeFilter;
}

VtkVisPipelineItem* VtkVisPipelineItem::child(int row) const
{
    TreeItem* treeItem = TreeItem::child(row);
    if (treeItem)
    {
        return dynamic_cast<VtkVisPipelineItem*>(treeItem);
    }

    return nullptr;
}

QVariant VtkVisPipelineItem::data(int column) const
{
    if (column == 1)
    {
        return isVisible();
    }

    return TreeItem::data(column);
}

bool VtkVisPipelineItem::setData(int column, const QVariant& value)
{
    if (column == 1)
    {
        setVisible(value.toBool());
        return true;
    }

    return TreeItem::setData(column, value);
}
bool VtkVisPipelineItem::isVisible() const
{
    return static_cast<bool>(_actor->GetVisibility());
}

void VtkVisPipelineItem::setVisible(bool visible)
{
    _actor->SetVisibility(static_cast<int>(visible));
    _actor->Modified();
    _renderer->Render();
}

int VtkVisPipelineItem::writeToFile(const std::string& filename) const
{
    if (!filename.empty())
    {
        return callVTKWriter(this->algorithm(), filename);
    }
    return 0;
}

int VtkVisPipelineItem::callVTKWriter(vtkAlgorithm* algorithm,
                                      const std::string& filename) const
{
    // needs to be implemented in derived classes!
    (void)algorithm;
    (void)filename;
    return 0;
}

vtkProp3D* VtkVisPipelineItem::actor() const
{
    return _actor;
}

void VtkVisPipelineItem::setScale(double x, double y, double z) const
{
    (void)x;
    (void)y, (void)z;
}

void VtkVisPipelineItem::setTranslation(double x, double y, double z) const
{
    (void)x;
    (void)y, (void)z;
}

void VtkVisPipelineItem::setScaleOnChildren(double x, double y, double z) const
{
    for (int i = 0; i < this->childCount(); ++i)
    {
        VtkVisPipelineItem* child = this->child(i);
        child->setScale(x, y, z);
    }
}

void VtkVisPipelineItem::setBackfaceCulling(bool enable) const
{
    // Reimplemented in subclass
    (void)enable;
}

void VtkVisPipelineItem::setBackfaceCullingOnChildren(bool enable) const
{
    for (int i = 0; i < this->childCount(); ++i)
    {
        VtkVisPipelineItem* child = this->child(i);
        child->setBackfaceCulling(static_cast<int>(enable));
        child->setBackfaceCullingOnChildren(static_cast<int>(enable));
    }
}

QStringList VtkVisPipelineItem::getScalarArrayNames() const
{
    this->algorithm()->Update();
    vtkDataSet* dataSet =
        vtkDataSet::SafeDownCast(this->algorithm()->GetOutputDataObject(0));
    QStringList dataSetAttributesList;
    if (dataSet)
    {
        vtkPointData* pointData = dataSet->GetPointData();
        // std::cout << "  #point data arrays: " <<
        // pointData->GetNumberOfArrays() << std::endl;
        for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
        {
            // std::cout << "    Name: " << pointData->GetArrayName(i) <<
            // std::endl;
            dataSetAttributesList.push_back(QString("P-") +
                                            pointData->GetArrayName(i));
        }

        vtkCellData* cellData = dataSet->GetCellData();
        // std::cout << "  #cell data arrays: " << cellData->GetNumberOfArrays()
        // << std::endl;
        for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
        {
            // std::cout << "    Name: " << cellData->GetArrayName(i) <<
            // std::endl;
            dataSetAttributesList.push_back(QString("C-") +
                                            cellData->GetArrayName(i));
        }
    }
    return dataSetAttributesList;
}
