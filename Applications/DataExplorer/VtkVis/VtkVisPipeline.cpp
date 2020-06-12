/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-17
 * \brief  Implementation of the VtkVisPipeline class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisPipeline.h"

#include <vtkAlgorithm.h>
#include <vtkCamera.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkImageActor.h>
#include <vtkImageReader2.h>
#include <vtkLight.h>
#include <vtkPointSet.h>
#include <vtkProp3D.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTransformFilter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>

#include <QColor>
#include <QFileInfo>
#include <QSettings>
#include <QString>
#include <QTime>

#include "BaseLib/Logging.h"

#include "MathLib/InterpolationAlgorithms/LinearIntervalInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"

#include "GeoTreeModel.h"
#include "MeshItem.h"
#include "MeshModel.h"
#include "StationTreeModel.h"
#include "TreeModel.h"
#include "VtkAlgorithmProperties.h"
#include "VtkCompositeElementSelectionFilter.h"
#include "VtkCompositeGeoObjectFilter.h"
#include "VtkCompositeNodeSelectionFilter.h"
#include "VtkFilterFactory.h"
#include "VtkVisImageItem.h"
#include "VtkVisPipelineItem.h"
#include "VtkVisPointSetItem.h"

VtkVisPipeline::VtkVisPipeline(vtkRenderer* renderer, QObject* parent /*= 0*/)
    : TreeModel(parent), renderer_(renderer)
{
    QList<QVariant> rootData;
    rootData << "Object name" << "Visible";
    delete rootItem_;
    rootItem_ = new TreeItem(rootData, nullptr);

    QSettings settings;
    QVariant backgroundColorVariant = settings.value("VtkBackgroundColor");
    if (backgroundColorVariant != QVariant())
    {
        this->setBGColor(backgroundColorVariant.value<QColor>());
    }

    resetCameraOnAddOrRemove_ = settings.value("resetViewOnLoad", true).toBool();
}

bool VtkVisPipeline::setData( const QModelIndex &index, const QVariant &value,
                              int role /* = Qt::EditRole */ )
{
    emit vtkVisPipelineChanged();

    return TreeModel::setData(index, value, role);
}

void VtkVisPipeline::addLight(const GeoLib::Point &pos)
{
    double lightPos[3];
    for (auto& light : lights_)
    {
        light->GetPosition(lightPos);
        if (pos[0] == lightPos[0] && pos[1] == lightPos[1] &&
            pos[2] == lightPos[2])
        {
            return;
        }
    }
    vtkLight* l = vtkLight::New();
    l->SetPosition(pos[0], pos[1], pos[2]);
    renderer_->AddLight(l);
    lights_.push_back(l);
}

vtkLight* VtkVisPipeline::getLight(const GeoLib::Point &pos) const
{
    double lightPos[3];
    for (auto light : lights_)
    {
        light->GetPosition(lightPos);
        if (pos[0] == lightPos[0] && pos[1] == lightPos[1] &&
            pos[2] == lightPos[2])
        {
            return light;
        }
    }
    return nullptr;
}

void VtkVisPipeline::removeLight(const GeoLib::Point &pos)
{
    double lightPos[3];
    for (auto it = lights_.begin(); it != lights_.end(); ++it)
    {
        (*it)->GetPosition(lightPos);
        if (pos[0] == lightPos[0] && pos[1] == lightPos[1] && pos[2] == lightPos[2])
        {
            renderer_->RemoveLight(*it);
            (*it)->Delete();
            lights_.erase(it);
            return;
        }
    }
}

QColor VtkVisPipeline::getBGColor() const
{
    double* color = renderer_->GetBackground();
    QColor c(static_cast<int>(color[0] * 255),
             static_cast<int>(color[1] * 255),
             static_cast<int>(color[2] * 255));
    return c;
}

void VtkVisPipeline::setBGColor(const QColor &color)
{
    QSettings settings;
    settings.setValue("VtkBackgroundColor", color);
    renderer_->SetBackground(color.redF(), color.greenF(), color.blueF());
}

QModelIndex VtkVisPipeline::getIndex( vtkProp3D* actor )
{
    return actorMap_.value(actor, QModelIndex());
}

Qt::ItemFlags VtkVisPipeline::flags( const QModelIndex &index ) const
{
    Qt::ItemFlags defaultFlags = Qt::ItemIsEnabled | Qt::ItemIsSelectable;

    if (!index.isValid())
    {
        return Qt::ItemIsEnabled;
    }

    //if (index.column() == 1)
    //    defaultFlags |= Qt::ItemIsEditable;

    return defaultFlags;
}

void VtkVisPipeline::loadFromFile(QString filename)
{
#ifndef NDEBUG
    QTime myTimer;
    myTimer.start();
    INFO("VTK Read: {:s}.", filename.toStdString());
#endif

    if (filename.size() > 0)
    {
        vtkSmartPointer<vtkXMLDataReader> reader;
        if (filename.endsWith("vti"))
        {
            reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        }
        else if (filename.endsWith("vtr"))
        {
            reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
        }
        else if (filename.endsWith("vts"))
        {
            reader = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
        }
        else if (filename.endsWith("vtp"))
        {
            reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        }
        else if (filename.endsWith("vtu"))
        {
            reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        }
        else if (filename.endsWith("vtk"))
        {
            vtkGenericDataObjectReader* oldStyleReader =
                    vtkGenericDataObjectReader::New();
            oldStyleReader->SetFileName(filename.toStdString().c_str());
            oldStyleReader->ReadAllFieldsOn();
            oldStyleReader->ReadAllScalarsOn();
            oldStyleReader->Update();
            vtkDataSet* dataSet = vtkDataSet::SafeDownCast(oldStyleReader->GetOutput());
            if (dataSet)
            {
                this->listArrays(dataSet);
                addPipelineItem(oldStyleReader);
            }
            else
                ERR("VtkVisPipeline::loadFromFile(): not a valid vtkDataSet.");
            return;
        }
        else
        {
            return;
        }

        reader->SetFileName(filename.toStdString().c_str());
        // TODO: insert ReadAllScalarsOn()-equivalent for xml-file-reader here, otherwise arrays are not available in GUI!
        reader->Update();
        //std::cout << "#cell scalars: " << reader->GetNumberOfCellArrays() << std::endl;
        //std::cout << "#point scalars: " << reader->GetNumberOfPointArrays() << std::endl;

        vtkSmartPointer<vtkDataSet> dataSet = reader->GetOutputAsDataSet();
        if (dataSet)
        {
            this->listArrays(dataSet);
            addPipelineItem(reader);
        }
        else
            ERR("VtkVisPipeline::loadFromFile(): not a valid vtkDataSet.");
    }

#ifndef NDEBUG
    INFO("{:d} ms", myTimer.elapsed());
#endif
}

void VtkVisPipeline::setGlobalSuperelevation(double factor) const
{
    // iterate over all source items
    for (int i = 0; i < rootItem_->childCount(); ++i)
    {
        auto* item = static_cast<VtkVisPipelineItem*>(rootItem_->child(i));
        item->setScale(1.0, 1.0, factor);

        // recursively set on all child items
        item->setScaleOnChildren(1.0, 1.0, 1.0);
    }

    emit vtkVisPipelineChanged();
}

void VtkVisPipeline::setGlobalBackfaceCulling(bool enable) const
{
    // iterate over all source items
    for (int i = 0; i < rootItem_->childCount(); ++i)
    {
        auto* item = static_cast<VtkVisPipelineItem*>(rootItem_->child(i));
        item->setBackfaceCulling(enable);

        // recursively set on all child items
        item->setBackfaceCullingOnChildren(enable);
    }

    emit vtkVisPipelineChanged();
}

void VtkVisPipeline::addPipelineItem(GeoTreeModel* model,
                                     const std::string &name,
                                     GeoLib::GEOTYPE type)
{
    addPipelineItem(model->vtkSource(name, type));
}

void VtkVisPipeline::addPipelineItem(StationTreeModel* model, const std::string &name)
{
    addPipelineItem(model->vtkSource(name));
}

void VtkVisPipeline::addPipelineItem(MeshModel* model, const QModelIndex &idx)
{
    addPipelineItem(static_cast<MeshItem*>(model->getItem(idx))->vtkSource());
}

QModelIndex VtkVisPipeline::addPipelineItem(VtkVisPipelineItem* item, const QModelIndex &parent)
{
    beginResetModel();

    item->Initialize(renderer_);
    TreeItem* parentItem = item->parentItem();
    parentItem->appendChild(item);

    if (!parent.isValid()) // Set global superelevation on source objects
    {
        QSettings settings;
        if (dynamic_cast<vtkImageAlgorithm*>(item->algorithm()) == nullptr)
        {  // if not an image
            item->setScale(
                1.0, 1.0,
                settings.value("globalSuperelevation", 1.0).toDouble());
        }
    }

    int parentChildCount = parentItem->childCount();
    QModelIndex newIndex = index(parentChildCount - 1, 0, parent);

    if (resetCameraOnAddOrRemove_ || rootItem_->childCount() == 1)
    {
        renderer_->ResetCamera(renderer_->ComputeVisiblePropBounds());
    }
    actorMap_.insert(item->actor(), newIndex);

    // Do not interpolate images
    if (dynamic_cast<vtkImageAlgorithm*>(item->algorithm()))
    {
        static_cast<vtkImageActor*>(item->actor())->InterpolateOff();
    }

    endResetModel();
    emit vtkVisPipelineChanged();
    emit itemSelected(newIndex);

    return newIndex;
}

QModelIndex VtkVisPipeline::addPipelineItem( vtkAlgorithm* source, QModelIndex parent /* = QModelindex() */)
{
    std::string itemName;

    if (!parent.isValid()) // if source object
    {
        auto* old_reader = dynamic_cast<vtkGenericDataObjectReader*>(source);
        auto* new_reader = dynamic_cast<vtkXMLReader*>(source);
        auto* image_reader = dynamic_cast<vtkImageReader2*>(source);
        auto* props = dynamic_cast<VtkAlgorithmProperties*>(source);
        auto* meshSource = dynamic_cast<MeshLib::VtkMappedMeshSource*>(source);
        if (old_reader)
        {
            itemName = old_reader->GetFileName();
        }
        else if (new_reader)
        {
            itemName = new_reader->GetFileName();
        }
        else if (image_reader)
        {
            itemName = image_reader->GetFileName();
        }
        else if (props)
        {
            itemName = props->GetName().toStdString();
        }
        else if (meshSource)
        {
            itemName = meshSource->GetMesh()->getName();
        }
    }

    if (itemName.length() == 0)
    {
        itemName = source->GetClassName();
    }

    QList<QVariant> itemData;
    itemData << QString::fromStdString(itemName) << true;

    VtkVisPipelineItem* item(nullptr);
    if (dynamic_cast<vtkImageAlgorithm*>(source))
    {
        item = new VtkVisImageItem(source, getItem(parent), itemData);
    }
    else
    {
        item = new VtkVisPointSetItem(source, getItem(parent), itemData);
    }
    return this->addPipelineItem(item, parent);
}

void VtkVisPipeline::removeSourceItem(GeoTreeModel* model,
                                      const std::string &name,
                                      GeoLib::GEOTYPE type)
{
    for (int i = 0; i < rootItem_->childCount(); i++)
    {
        VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(getItem(index(i, 0)));
        if (item->algorithm() == model->vtkSource(name, type))
        {
            removePipelineItem(index(i, 0));
            return;
        }
    }
}

void VtkVisPipeline::removeSourceItem(StationTreeModel* model, const std::string &name)
{
    for (int i = 0; i < rootItem_->childCount(); i++)
    {
        VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(getItem(index(i, 0)));
        if (item->algorithm() == model->vtkSource(name))
        {
            removePipelineItem(index(i, 0));
            return;
        }
    }
}

void VtkVisPipeline::removeSourceItem(MeshModel* model, const QModelIndex &idx)
{
    auto* sItem = static_cast<MeshItem*>(model->getItem(idx));

    for (int i = 0; i < rootItem_->childCount(); i++)
    {
        VtkVisPipelineItem* item =
            static_cast<VtkVisPipelineItem*>(getItem(index(i, 0)));
        if (item->algorithm() == sItem->vtkSource())
        {
            removePipelineItem(index(i, 0));
            return;
        }
    }
}

void VtkVisPipeline::removePipelineItem( QModelIndex index )
{
    if (!index.isValid())
    {
        return;
    }

    QMap<vtkProp3D*, QModelIndex>::iterator it = actorMap_.begin();
    while (it != actorMap_.end())
    {
        QModelIndex itIndex = it.value();
        if (itIndex == index)
        {
            actorMap_.erase(it);
            break;
        }
        ++it;
    }

    //TreeItem* item = getItem(index);
    removeRows(index.row(), 1, index.parent());

    if (resetCameraOnAddOrRemove_)
    {
        renderer_->ResetCamera(renderer_->ComputeVisiblePropBounds());
    }
    emit vtkVisPipelineChanged();
}

void VtkVisPipeline::listArrays(vtkDataSet* dataSet)
{
    if (dataSet)
    {
        vtkPointData* pointData = dataSet->GetPointData();
        INFO("  #point data arrays: {:d}", pointData->GetNumberOfArrays());
        for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
            INFO("    Name: {:s}", pointData->GetArrayName(i));

        vtkCellData* cellData = dataSet->GetCellData();
        INFO("  #cell data arrays: {:d}", cellData->GetNumberOfArrays());
        for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
            INFO("    Name: {:s}", cellData->GetArrayName(i));
    }
    else
        ERR("VtkVisPipeline::listArrays(): not a valid vtkDataSet.");
}

void VtkVisPipeline::showMeshElementQuality(
    MeshLib::VtkMappedMeshSource* source,
    MeshLib::MeshQualityType t, std::vector<double> const& quality)
{
    if (!source || quality.empty())
    {
        return;
    }

    int const nSources = this->rootItem_->childCount();
    for (int i = 0; i < nSources; i++)
    {
        auto* parentItem =
            static_cast<VtkVisPipelineItem*>(rootItem_->child(i));
        if (parentItem->algorithm() != source)
        {
            continue;
        }

        QList<QVariant> itemData;
        itemData << "MeshQuality: " +
                        QString::fromStdString(MeshQualityType2String(t))
                 << true;

        VtkCompositeFilter* filter = VtkFilterFactory::CreateCompositeFilter(
            "VtkCompositeElementSelectionFilter",
            parentItem->transformFilter());
        if (t == MeshLib::MeshQualityType::ELEMENTSIZE)
        {
            auto const range(
                std::minmax_element(quality.cbegin(), quality.cend()));
            static_cast<VtkCompositeElementSelectionFilter*>(filter)->setRange(
                *range.first, *range.second);
        }
        static_cast<VtkCompositeElementSelectionFilter*>(filter)->setSelectionArray("Selection", quality);
        VtkVisPointSetItem* item = new VtkVisPointSetItem(filter, parentItem, itemData);
        this->addPipelineItem(item, this->createIndex(i, 0, item));
        break;
    }
}

void VtkVisPipeline::highlightGeoObject(const vtkPolyDataAlgorithm* source, int index)
{
    this->removeHighlightedGeoObject();
    int nSources = this->rootItem_->childCount();
    for (int i = 0; i < nSources; i++)
    {
        auto* parentItem =
            static_cast<VtkVisPipelineItem*>(rootItem_->child(i));
        if (parentItem->algorithm() == source)
        {
            QList<QVariant> itemData;
            itemData << "Selected GeoObject" << true;

            VtkCompositeFilter* filter =
                VtkFilterFactory::CreateCompositeFilter(
                    "VtkCompositeGeoObjectFilter",
                    parentItem->transformFilter());
            static_cast<VtkCompositeGeoObjectFilter*>(filter)->SetIndex(index);
            VtkVisPointSetItem* item =
                new VtkVisPointSetItem(filter, parentItem, itemData);
            QModelIndex parent_index =
                static_cast<TreeModel*>(this)->index(i, 0, QModelIndex());
            highlighted_geo_index_ = this->addPipelineItem(item, parent_index);
            break;
        }
    }
}

void VtkVisPipeline::removeHighlightedGeoObject()
{
    if (highlighted_geo_index_ != QModelIndex())
    {
        this->removePipelineItem(highlighted_geo_index_);
        highlighted_geo_index_ = QModelIndex();
    }
}

void VtkVisPipeline::highlightMeshComponent(vtkUnstructuredGridAlgorithm const*const source, unsigned index, bool is_element)
{
    int nSources = this->rootItem_->childCount();
    for (int i = 0; i < nSources; i++)
    {
        auto* parentItem =
            static_cast<VtkVisPipelineItem*>(rootItem_->child(i));
        if (parentItem->algorithm() == source)
        {
            QList<QVariant> itemData;
            itemData << "Selected Mesh Component" << true;
            QList<QVariant> selected_index;
            selected_index << index << index;

            VtkCompositeFilter* filter(nullptr);
            if (is_element)
            {
                filter = VtkFilterFactory::CreateCompositeFilter(
                    "VtkCompositeElementSelectionFilter",
                    parentItem->transformFilter());
                static_cast<VtkCompositeElementSelectionFilter*>(filter)
                    ->setSelectionArray("vtkIdFilter_Ids");
                static_cast<VtkCompositeElementSelectionFilter*>(filter)
                    ->SetUserVectorProperty("Threshold Between",
                                            selected_index);
            }
            else
            {
                filter = VtkFilterFactory::CreateCompositeFilter(
                    "VtkCompositeNodeSelectionFilter",
                    parentItem->transformFilter());
                std::vector<unsigned> indeces(1);
                indeces[0] = index;
                static_cast<VtkCompositeNodeSelectionFilter*>(filter)
                    ->setSelectionArray(indeces);
            }
            VtkVisPointSetItem* item =
                new VtkVisPointSetItem(filter, parentItem, itemData);
            QModelIndex parent_index =
                static_cast<TreeModel*>(this)->index(i, 0, QModelIndex());
            highlighted_mesh_component_ =
                this->addPipelineItem(item, parent_index);
            break;
        }
    }
}

void VtkVisPipeline::removeHighlightedMeshComponent()
{
    if (highlighted_mesh_component_ != QModelIndex())
    {
        this->removePipelineItem(highlighted_mesh_component_);
        highlighted_mesh_component_ = QModelIndex();
    }
}
