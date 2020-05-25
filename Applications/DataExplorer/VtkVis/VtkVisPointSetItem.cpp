/**
 * \file
 * \author Karsten Rink
 * \date   2011-09-29
 * \brief  Implementation of the VtkVisPointSetItem class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisPointSetItem.h"

#include <limits>

#include <QObject>
#include <QRegExp>
#include <QSettings>
#include <QStringList>

#include "BaseLib/Logging.h"

#include "BaseLib/FileTools.h"

#include "VtkAlgorithmProperties.h"
#include "VtkCompositeFilter.h"
#include "VtkCompositeContourFilter.h"
#include "VtkCompositeThresholdFilter.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"

#include "QVtkDataSetMapper.h"
#include <vtkActor.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>
#include <vtkImageAlgorithm.h>
#include <vtkPointData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>

// export test
#include <vtkPolyDataAlgorithm.h>
#include <vtkTriangleFilter.h>
#include <vtkTubeFilter.h>
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkDataSetAttributes.h>

VtkVisPointSetItem::VtkVisPointSetItem(
        vtkAlgorithm* algorithm, TreeItem* parentItem,
        const QList<QVariant> data /*= QList<QVariant>()*/)
    : VtkVisPipelineItem(algorithm, parentItem, data), mapper_(nullptr),
    transformFilter_(nullptr), onPointData_(true), activeArrayName_("")
{
    auto* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
    if (parentItem->parentItem())
    {
        // special case if parent is image but child is not (e.g.
        // Image2BarChartFilter)
        if (dynamic_cast<vtkImageAlgorithm*>(visParentItem->algorithm()))
        {
            algorithm_->SetInputConnection(
                visParentItem->algorithm()->GetOutputPort());
        }
        else
        {
            auto* pointSetItem = dynamic_cast<VtkVisPointSetItem*>(parentItem);
            if (pointSetItem)
            {
                algorithm_->SetInputConnection(
                    pointSetItem->transformFilter()->GetOutputPort());
            }
        }
    }
}

VtkVisPointSetItem::VtkVisPointSetItem(
        VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
        const QList<QVariant> data /*= QList<QVariant>()*/)
    : VtkVisPipelineItem(compositeFilter, parentItem, data), mapper_(nullptr),
    transformFilter_(nullptr), onPointData_(true), activeArrayName_("")
{
}

VtkVisPointSetItem::~VtkVisPointSetItem()
{
    transformFilter_->Delete();
    mapper_->Delete();
}
QString VtkVisPointSetItem::GetActiveAttribute() const
{
    return vtkProps_->GetActiveAttribute();
}

void VtkVisPointSetItem::Initialize(vtkRenderer* renderer)
{
    // TODO vtkTransformFilter creates a new copy of the point coordinates which
    // conflicts with VtkMappedMeshSource. Find a workaround!
    transformFilter_ = vtkTransformFilter::New();
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Identity();
    transformFilter_->SetTransform(transform);

    transformFilter_->SetInputConnection(algorithm_->GetOutputPort());
    transformFilter_->Update();

    renderer_ = renderer;
    mapper_ = QVtkDataSetMapper::New();
    mapper_->InterpolateScalarsBeforeMappingOff();
    mapper_->SetColorModeToMapScalars();

    mapper_->SetInputConnection(transformFilter_->GetOutputPort());
    actor_ = vtkActor::New();
    static_cast<vtkActor*>(actor_)->SetMapper(mapper_);
    renderer_->AddActor(actor_);

    // Determine the right pre-set properties
    // Order is: algorithm_, compositeFilter_, create a new one with props copied from parent
    auto* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(algorithm_);
    if (!vtkProps)
    {
        vtkProps = dynamic_cast<VtkAlgorithmProperties*>(compositeFilter_);

        // Copy properties from parent or create a new VtkAlgorithmProperties
        if (!vtkProps)
        {
            auto* parentItem =
                dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
            while (parentItem)
            {
                VtkAlgorithmProperties* parentProps = nullptr;
                if (dynamic_cast<VtkVisPointSetItem*>(parentItem))
                {
                    parentProps = dynamic_cast<VtkVisPointSetItem*>(parentItem)
                                      ->getVtkProperties();
                }
                if (parentProps)
                {
                    vtkProps =
                        new VtkAlgorithmProperties();  // TODO memory leak?
                    vtkProps->SetScalarVisibility(
                        parentProps->GetScalarVisibility());
                    vtkProps->SetTexture(parentProps->GetTexture());
                    vtkProps->SetActiveAttribute(
                        parentProps->GetActiveAttribute());
                    parentItem = nullptr;
                }
                else
                {
                    parentItem = dynamic_cast<VtkVisPipelineItem*>(
                        parentItem->parentItem());
                }
            }

            // Has no parents
            if (!vtkProps)
            {
                vtkProps = new VtkAlgorithmProperties();  // TODO memory leak?
            }
        }
    }
    vtkProps_ = vtkProps;

    if (vtkProps->GetActiveAttribute().length() == 0)
    {
        // Get first scalar and set it to active
        QStringList arrayNames = this->getScalarArrayNames();
        if (arrayNames.length() > 0)
        {
            vtkProps->SetActiveAttribute(arrayNames[0]);
        }
        else
        {
            vtkProps->SetActiveAttribute("Solid Color");
        }
    }
    this->setVtkProperties(vtkProps);
    this->SetActiveAttribute(vtkProps->GetActiveAttribute());


    // Set global backface culling
    QSettings settings;
    bool backfaceCulling = settings.value("globalCullBackfaces", 0).toBool();
    this->setBackfaceCulling(backfaceCulling);

    // Set the correct threshold range
    if (dynamic_cast<VtkCompositeThresholdFilter*>(this->compositeFilter_) )
    {
        double range[2];
        this->GetRangeForActiveAttribute(range);
        QList<QVariant> thresholdRangeList;
        thresholdRangeList.push_back(range[0]);
        thresholdRangeList.push_back(range[1]);
        dynamic_cast<VtkCompositeFilter*>(this->compositeFilter_)
            ->SetUserVectorProperty("Range", thresholdRangeList);
    }

    // Show edges on meshes
    if (dynamic_cast<MeshLib::VtkMappedMeshSource*>(this->algorithm_))
    {
        vtkProps_->GetProperties()->SetEdgeVisibility(1);
    }
}

void VtkVisPointSetItem::SetScalarVisibility( bool on )
{
    mapper_->SetScalarVisibility(on);
}

void VtkVisPointSetItem::setVtkProperties(VtkAlgorithmProperties* vtkProps)
{
    QObject::connect(vtkProps, SIGNAL(ScalarVisibilityChanged(bool)),
                     mapper_, SLOT(SetScalarVisibility(bool)));

    auto* actor = dynamic_cast<vtkActor*>(actor_);
    if (actor)
    {
        if (vtkProps->GetTexture() != nullptr)
        {
            vtkProps->SetScalarVisibility(false);
            actor->GetProperty()->SetColor(1,1,1); // don't colorise textures
            actor->SetTexture(vtkProps->GetTexture());
        }
        else
        {
            vtkSmartPointer<vtkProperty> itemProperty = vtkProps->GetProperties();
            actor->SetProperty(itemProperty);
        }

        if (!vtkProps->GetScalarVisibility())
        {
            vtkProps->SetScalarVisibility(false);
        }
    }
}

int VtkVisPointSetItem::callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const
{
    std::string file_name_cpy(filename);
    auto* algPD = dynamic_cast<vtkPolyDataAlgorithm*>(algorithm);
    if (algPD)
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> pdWriter =
                vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        pdWriter->SetInputData(algPD->GetOutputDataObject(0));
        if (BaseLib::getFileExtension(filename) != ".vtp")
        {
            file_name_cpy.append(".vtp");
        }
        pdWriter->SetFileName(file_name_cpy.c_str());
        return pdWriter->Write();
    }

    auto* algUG = dynamic_cast<vtkUnstructuredGridAlgorithm*>(algorithm);
    if (algUG)
    {
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> ugWriter =
                vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        ugWriter->SetInputData(algUG->GetOutputDataObject(0));
        if (BaseLib::getFileExtension(filename) != ".vtu")
        {
            file_name_cpy.append(".vtu");
        }
        ugWriter->SetFileName(file_name_cpy.c_str());
        return ugWriter->Write();
    }

    WARN("VtkVisPipelineItem::writeToFile(): Unknown data type.");
    return 0;
}

void VtkVisPointSetItem::SetActiveAttribute( const QString& name )
{
    // Get type by identifier
    if (name.contains(QRegExp("^P-")))
    {
        onPointData_ = true;
    }
    else if (name.contains(QRegExp("^C-")))
    {
        onPointData_ = false;
    }
    else if (name.contains("Solid Color"))
    {
        vtkProps_->SetActiveAttribute("Solid Color");
        mapper_->ScalarVisibilityOff();
        return;
    }
    else
    {
        return;
    }

    // Remove type identifier
    activeArrayName_ = QString(name).remove(0, 2).toStdString();

    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->algorithm_->GetOutputDataObject(0));
    if (!dataSet)
    {
        return;
    }

    double range[2];
    GetRangeForActiveAttribute(range);
    if (onPointData_)
    {
        vtkPointData* pointData = dataSet->GetPointData();
        if(pointData)
        {
            algorithm_->SetInputArrayToProcess(0, 0, 0,
                vtkDataObject::FIELD_ASSOCIATION_POINTS, activeArrayName_.c_str());
            mapper_->SetScalarModeToUsePointFieldData();
        }
    }
    else
    {
        vtkCellData* cellData = dataSet->GetCellData();
        if(cellData)
        {
            algorithm_->SetInputArrayToProcess(0, 0, 0,
                vtkDataObject::FIELD_ASSOCIATION_CELLS, activeArrayName_.c_str());
            mapper_->SetScalarModeToUseCellFieldData();
        }
    }

    vtkProps_->SetActiveAttribute(name);

    mapper_->ScalarVisibilityOn();
    mapper_->UseLookupTableScalarRangeOn();

    auto* mapper = dynamic_cast<QVtkDataSetMapper*>(mapper_);
    if (mapper)
    {
        // Create a default color table when there is no lookup table for this attribute
        vtkLookupTable* lut = vtkProps_->GetLookupTable(name);
        if (lut == nullptr)
        {
            //std::cout << "Creating new lookup table for: " << name.toStdString() << std::endl;
            lut = vtkLookupTable::New(); // is not a memory leak, gets deleted in VtkAlgorithmProperties
            lut->SetTableRange(range);
            vtkProps_->SetLookUpTable(name, lut);
        }
        mapper_->SetLookupTable(lut);
    }
    mapper_->SelectColorArray( activeArrayName_.c_str());
}

bool VtkVisPointSetItem::activeAttributeExists(vtkDataSetAttributes* data, std::string& name)
{
    bool arrayFound = false;
    for (int i = 0; i < data->GetNumberOfArrays() && !arrayFound; i++)
    {
        std::string arrayName = data->GetArrayName(i);
        if (arrayName == name)
        {
            arrayFound = true;
        }
    }
    if (arrayFound)
    {
        // TODO Necessary? Currently this function is not called
        data->SetActiveAttribute(name.c_str(), vtkDataSetAttributes::SCALARS);
        return true;
    }

    return false;
}

void VtkVisPointSetItem::setScale(double x, double y, double z) const
{
    if (this->transformFilter())
    {
        auto* transform =
            static_cast<vtkTransform*>(this->transformFilter_->GetTransform());
        double* trans = transform->GetPosition();
        transform->Identity();
        transform->Scale(x, y, z);
        transform->Translate(trans[0] / x, trans[1] / y, trans[2] / z);
        this->transformFilter()->Modified();
    }
}

void VtkVisPointSetItem::setTranslation(double x, double y, double z) const
{
    if (this->transformFilter())
    {
        auto* transform =
            static_cast<vtkTransform*>(this->transformFilter_->GetTransform());
        double* scale = transform->GetScale();
        transform->Identity();
        transform->Scale(scale);
        transform->Translate(x, y, z);
        this->transformFilter()->Modified();
    }
}

vtkAlgorithm* VtkVisPointSetItem::transformFilter() const
{
    return transformFilter_;
}

void VtkVisPointSetItem::setBackfaceCulling(bool enable) const
{
    static_cast<vtkActor*>(this->actor_)
        ->GetProperty()
        ->SetBackfaceCulling(static_cast<int>(enable));
}

void VtkVisPointSetItem::GetRangeForActiveAttribute(double range[2]) const
{
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->algorithm_->GetOutputDataObject(0));
    if (dataSet && activeArrayName_.length() > 0)
    {
        if (onPointData_)
        {
            vtkPointData* pointData = dataSet->GetPointData();
            if (pointData)
            {
                pointData->GetArray(activeArrayName_.c_str())->GetRange(range);
            }
        }
        else
        {
            vtkCellData* cellData = dataSet->GetCellData();
                cellData->GetArray(activeArrayName_.c_str())->GetRange(range);
        }
    }
}
