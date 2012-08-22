/**
 * \file VtkVisPointSetItem.cpp
 * 2011/09/29 KR Initial implementation
 */

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include "VtkVisPointSetItem.h"
#include "VtkCompositeFilter.h"
#include "VtkCompositeContourFilter.h"
#include "VtkCompositeThresholdFilter.h"

#include <limits>

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

#include <QObject>
#include <QRegExp>
#include <QSettings>
#include <QStringList>

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
	: VtkVisPipelineItem(algorithm, parentItem, data), _mapper(NULL),
	_transformFilter(NULL), _onPointData(true), _activeArrayName("")
{
	VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
	if (parentItem->parentItem())
	{
		// special case if parent is image but child is not (e.g. Image2BarChartFilter)
		if (dynamic_cast<vtkImageAlgorithm*>(visParentItem->algorithm()))
			_algorithm->SetInputConnection(visParentItem->algorithm()->GetOutputPort());
		else
		{
			VtkVisPointSetItem* pointSetItem =
			        dynamic_cast<VtkVisPointSetItem*>(parentItem);
			if (pointSetItem)
				_algorithm->SetInputConnection(
				        pointSetItem->transformFilter()->GetOutputPort());
		}
	}
}

VtkVisPointSetItem::VtkVisPointSetItem(
        VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
        const QList<QVariant> data /*= QList<QVariant>()*/)
	: VtkVisPipelineItem(compositeFilter, parentItem, data), _mapper(NULL),
	_transformFilter(NULL), _onPointData(true), _activeArrayName("")
{
}

VtkVisPointSetItem::~VtkVisPointSetItem()
{
	_transformFilter->Delete();
	_mapper->Delete();
}
const QString VtkVisPointSetItem::GetActiveAttribute() const
{
	return _vtkProps->GetActiveAttribute();
}

void VtkVisPointSetItem::Initialize(vtkRenderer* renderer)
{
	_transformFilter = vtkTransformFilter::New();
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Identity();
	_transformFilter->SetTransform(transform);

	_transformFilter->SetInputConnection(_algorithm->GetOutputPort());
	_transformFilter->Update();

	_renderer = renderer;
	_mapper = QVtkDataSetMapper::New();
	_mapper->InterpolateScalarsBeforeMappingOff();
	_mapper->SetColorModeToMapScalars();

	_mapper->SetInputConnection(_transformFilter->GetOutputPort());
	_actor = vtkActor::New();
	static_cast<vtkActor*>(_actor)->SetMapper(_mapper);
	_renderer->AddActor(_actor);

	// Determine the right pre-set properties
	// Order is: _algorithm, _compositeFilter, create a new one with props copied from parent
	VtkAlgorithmProperties* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(_algorithm);
	if (!vtkProps)
	{
		vtkProps = dynamic_cast<VtkAlgorithmProperties*>(_compositeFilter);

		// Copy properties from parent or create a new VtkAlgorithmProperties
		if (!vtkProps)
		{
			VtkVisPipelineItem* parentItem = dynamic_cast<VtkVisPipelineItem*>(this->parentItem());
			while (parentItem)
			{
				VtkAlgorithmProperties* parentProps = NULL;
				if(dynamic_cast<VtkVisPointSetItem*>(parentItem))
					parentProps = dynamic_cast<VtkVisPointSetItem*>(parentItem)->getVtkProperties();
				if (parentProps)
				{
					vtkProps = new VtkAlgorithmProperties(); // TODO memory leak?
					vtkProps->SetScalarVisibility(parentProps->GetScalarVisibility());
					vtkProps->SetTexture(parentProps->GetTexture());
					vtkProps->SetActiveAttribute(parentProps->GetActiveAttribute());
					parentItem = NULL;
				}
				else
					parentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem->parentItem());
			}

			// Has no parents
			if (!vtkProps)
				vtkProps = new VtkAlgorithmProperties(); // TODO memory leak?
		}
	}
	_vtkProps = vtkProps;

	if (vtkProps->GetActiveAttribute().length() == 0)
	{
		// Get first scalar and set it to active
		QStringList arrayNames = this->getScalarArrayNames();
		if (arrayNames.length() > 0)
			vtkProps->SetActiveAttribute(arrayNames[0]);
		else
			vtkProps->SetActiveAttribute("Solid Color");
	}
	this->setVtkProperties(vtkProps);
	this->SetActiveAttribute(vtkProps->GetActiveAttribute());


	// Set global backface culling
	QSettings settings("UFZ, OpenGeoSys-5");
	bool backfaceCulling = settings.value("globalCullBackfaces", 0).toBool();
	this->setBackfaceCulling(backfaceCulling);

	// Set the correct threshold range
	if ( dynamic_cast<VtkCompositeThresholdFilter*>(this->_compositeFilter) )
	{
		double range[2];
		this->GetRangeForActiveAttribute(range);
		QList<QVariant> thresholdRangeList;
		thresholdRangeList.push_back(range[0]);
		thresholdRangeList.push_back(range[1]);
		dynamic_cast<VtkCompositeFilter*>(this->_compositeFilter)
			->SetUserVectorProperty("Range", thresholdRangeList);
	}
}

void VtkVisPointSetItem::SetScalarVisibility( bool on )
{
	_mapper->SetScalarVisibility(on);
}

void VtkVisPointSetItem::setVtkProperties(VtkAlgorithmProperties* vtkProps)
{
	QObject::connect(vtkProps, SIGNAL(ScalarVisibilityChanged(bool)),
	                 _mapper, SLOT(SetScalarVisibility(bool)));

	vtkActor* actor = dynamic_cast<vtkActor*>(_actor);
	if (actor)
	{
		if (vtkProps->GetTexture() != NULL)
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
			vtkProps->SetScalarVisibility(false);
	}
}

int VtkVisPointSetItem::callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const
{
	vtkPolyDataAlgorithm* algPD = dynamic_cast<vtkPolyDataAlgorithm*>(algorithm);
	vtkUnstructuredGridAlgorithm* algUG = dynamic_cast<vtkUnstructuredGridAlgorithm*>(algorithm);
	if (algPD)
	{
//		vtkGenericDataObjectWriter* pdWriter = vtkGenericDataObjectWriter::New();
		vtkSmartPointer<vtkXMLPolyDataWriter> pdWriter =
		        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		pdWriter->SetInput(algPD->GetOutputDataObject(0));
		//pdWriter->SetDataModeToAscii();
		//pdWriter->SetCompressorTypeToNone();
		std::string filenameWithExt = filename;
		filenameWithExt.append(".vtp");
		pdWriter->SetFileName(filenameWithExt.c_str());
		return pdWriter->Write();
	}
	else if (algUG)
	{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> ugWriter =
		        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		ugWriter->SetInput(algUG->GetOutputDataObject(0));
		//ugWriter->SetDataModeToAscii();
		//ugWriter->SetCompressorTypeToNone();
		std::string filenameWithExt = filename;
		filenameWithExt.append(".vtu");
		ugWriter->SetFileName(filenameWithExt.c_str());
		return ugWriter->Write();
	}
	std::cout << "VtkVisPipelineItem::writeToFile() - Unknown data type..." << std::endl;
	return 0;
}

void VtkVisPointSetItem::SetActiveAttribute( const QString& name )
{
	// Get type by identifier
	if (name.contains(QRegExp("^P-")))
		_onPointData = true;
	else if (name.contains(QRegExp("^C-")))
		_onPointData = false;
	else if (name.contains("Solid Color"))
	{
		_vtkProps->SetActiveAttribute("Solid Color");
		_mapper->ScalarVisibilityOff();
		return;
	}
	else
		return;

	// Remove type identifier
	_activeArrayName = QString(name).remove(0, 2).toStdString();
	const char* charName = _activeArrayName.c_str();

	double range[2];
	vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->_algorithm->GetOutputDataObject(0));
	if (dataSet)
	{
		if (_onPointData)
		{
			vtkPointData* pointData = dataSet->GetPointData();
			if(pointData)
			{
				if(activeAttributeExists(pointData, _activeArrayName))
				{
					_algorithm->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, charName);
					_mapper->SetScalarModeToUsePointData();
					pointData->GetArray(_activeArrayName.c_str())->GetRange(range);
				}
				else
				{
					_activeArrayName = "";
					_vtkProps->SetActiveAttribute("Solid Color");
					_mapper->ScalarVisibilityOff();
					return;
				}
			}
		}
		else
		{
			vtkCellData* cellData = dataSet->GetCellData();
			if(cellData)
			{
				if(activeAttributeExists(cellData, _activeArrayName))
				{
					_algorithm->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, charName);
					_mapper->SetScalarModeToUseCellData();
					cellData->GetArray(_activeArrayName.c_str())->GetRange(range);
				}
				else
				{
					_activeArrayName = "";
					_vtkProps->SetActiveAttribute("Solid Color");
					_mapper->ScalarVisibilityOff();
					return;
				}
			}
		}

		//std::cout << "Range for " << name.toStdString() << " :" << range[0] << " " << range[1] << std::endl;
		_vtkProps->SetActiveAttribute(name);

		QVtkDataSetMapper* mapper = dynamic_cast<QVtkDataSetMapper*>(_mapper);
		if (mapper)
		{
			// Create a default color table when there is no lookup table for this attribute
			vtkLookupTable* lut = _vtkProps->GetLookupTable(name);
			if (lut == NULL)
			{
				//std::cout << "Creating new lookup table for: " << name.toStdString() << std::endl;
				lut = vtkLookupTable::New(); // is not a memory leak, gets deleted in VtkAlgorithmProperties
				lut->SetTableRange(range);
				_vtkProps->SetLookUpTable(name, lut);
			}

			_mapper->SetLookupTable(lut);
			_mapper->UseLookupTableScalarRangeOn();
			//_mapper->SetScalarRange(range); // not necessary when UseLookupTableScalarRange is on
		}

		_mapper->ScalarVisibilityOn();
		_mapper->Update();
	}
}

bool VtkVisPointSetItem::activeAttributeExists(vtkDataSetAttributes* data, std::string& name)
{
	bool arrayFound = false;
	for (int i = 0; i < data->GetNumberOfArrays() && !arrayFound; i++)
	{
		std::string arrayName = data->GetArrayName(i);
		if(arrayName.compare(name) == 0)
			arrayFound = true;
	}
	if(arrayFound)
	{
		data->SetActiveAttribute(name.c_str(), vtkDataSetAttributes::SCALARS);
		return true;
	}
	else
		return false;
}

void VtkVisPointSetItem::setScale(double x, double y, double z) const
{
	if (this->transformFilter())
	{
		vtkTransform* transform =
			static_cast<vtkTransform*>(this->_transformFilter->GetTransform());
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
		vtkTransform* transform =
			static_cast<vtkTransform*>(this->_transformFilter->GetTransform());
		double* scale = transform->GetScale();
		transform->Identity();
		transform->Scale(scale);
		transform->Translate(x, y, z);
		this->transformFilter()->Modified();
	}
}

vtkAlgorithm* VtkVisPointSetItem::transformFilter() const
{
	return _transformFilter;
}

void VtkVisPointSetItem::setBackfaceCulling(bool enable) const
{
	static_cast<vtkActor*>(this->_actor)->GetProperty()->SetBackfaceCulling((int)enable);
}

void VtkVisPointSetItem::GetRangeForActiveAttribute(double range[2]) const
{
	vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->_algorithm->GetOutputDataObject(0));
	if (dataSet && _activeArrayName.length() > 0)
	{
		if (_onPointData)
		{
			vtkPointData* pointData = dataSet->GetPointData();
			if(pointData)
				pointData->GetArray(_activeArrayName.c_str())->GetRange(range);
		}
		else
		{
			vtkCellData* cellData = dataSet->GetCellData();
				cellData->GetArray(_activeArrayName.c_str())->GetRange(range);
		}
	}
}
