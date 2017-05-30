/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-17
 * \brief  Implementation of the VtkVisPipelineItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "FileTools.h"
#include "VtkAlgorithmProperties.h"
#include "VtkVisPipelineItem.h"
#include "VtkCompositeFilter.h"

#include "QVtkDataSetMapper.h"
#include <vtkActor.h>
#include <vtkAlgorithm.h>
#include <vtkDataSetMapper.h>
#include <vtkPointSet.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTextureMapToPlane.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkImageActor.h>

#include <QMessageBox>

#ifdef VTKFBXCONVERTER_FOUND
#include "ThirdParty/VtkFbxConverter/VtkFbxConverter.h"
#include "Common.h"
#include <fbxsdk.h>


extern FbxManager* lSdkManager;
extern FbxScene* lScene;
#endif

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
        _algorithm->SetInputConnection(
            visParentItem->algorithm()->GetOutputPort());
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

VtkVisPipelineItem* VtkVisPipelineItem::child( int row ) const
{
    TreeItem* treeItem = TreeItem::child(row);
    if (treeItem)
        return dynamic_cast<VtkVisPipelineItem*>(treeItem);

    return nullptr;
}

QVariant VtkVisPipelineItem::data( int column ) const
{
    if (column == 1)
        return isVisible();

    return TreeItem::data(column);
}

bool VtkVisPipelineItem::setData( int column, const QVariant &value )
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
    return (bool)_actor->GetVisibility();
}

void VtkVisPipelineItem::setVisible( bool visible )
{
    _actor->SetVisibility((int)visible);
    _actor->Modified();
    _renderer->Render();
}

int VtkVisPipelineItem::writeToFile(const std::string &filename) const
{
    if (!filename.empty())
    {
#ifdef VTKFBXCONVERTER_FOUND
        if (filename.substr(filename.size() - 4).find("fbx") != std::string::npos)
        {
            if(!dynamic_cast<vtkImageActor*>(_actor))
            {
                InitializeSdkObjects(lSdkManager, lScene);

                VtkFbxConverter fbxConverter(static_cast<vtkActor*>(_actor), lScene);
                fbxConverter.convert(BaseLib::extractBaseNameWithoutExtension(filename));
                FbxNode* node = fbxConverter.getNode();
                if(node)
                {
                    fbxConverter.addUserProperty("UseVertexColors", _vtkProps->GetScalarVisibility());
                    lScene->GetRootNode()->AddChild(node);
                    // Get the file format. Use either "FBX [6.0] binary (*.fbx)" or "FBX [6.0] ascii (*.fbx)"
                    int fbxFormat = lSdkManager->GetIOPluginRegistry()
                        ->FindWriterIDByDescription("FBX 6.0 binary (*.fbx)");
                    // Embed only works in "FBX 6.0 binary (*.fbx)"
                    const bool fbxEmbed = true;
                    SaveScene(lSdkManager, lScene, filename.c_str(), fbxFormat, fbxEmbed);
                    lScene->Clear();
                }
            }
            else
                QMessageBox::warning(
                    nullptr, "Conversion to FBX not possible",
                    "It is not possible to convert an vtkImageData based object \
                    to FBX. If you want to convert raster data import it via \" \
                    File / Import / Raster Files as PolyData\"!");
            return 0;
        }
#endif // VTKFBXCONVERTER_FOUND

        return callVTKWriter(this->algorithm(), filename);
    }
    return 0;
}

int VtkVisPipelineItem::callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const
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
        child->setBackfaceCulling((int)enable);
        child->setBackfaceCullingOnChildren((int)enable);
    }
}

QStringList VtkVisPipelineItem::getScalarArrayNames() const
{
    this->algorithm()->Update();
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->algorithm()->GetOutputDataObject(0));
    QStringList dataSetAttributesList;
    if (dataSet)
    {
        vtkPointData* pointData = dataSet->GetPointData();
        //std::cout << "  #point data arrays: " << pointData->GetNumberOfArrays() << std::endl;
        for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
        {
            // std::cout << "    Name: " << pointData->GetArrayName(i) << std::endl;
            dataSetAttributesList.push_back(QString("P-") + pointData->GetArrayName(i));
        }

        vtkCellData* cellData = dataSet->GetCellData();
        //std::cout << "  #cell data arrays: " << cellData->GetNumberOfArrays() << std::endl;
        for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
        {
            // std::cout << "    Name: " << cellData->GetArrayName(i) << std::endl;
            dataSetAttributesList.push_back(QString("C-") + cellData->GetArrayName(i));
        }
    }
    return dataSetAttributesList;
}
