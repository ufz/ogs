/**
 * \file VtkVisPointSetItem.cpp
 * 2011/09/29 KR Initial implementation
 */

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include "VtkVisPointSetItem.h"

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

#include <QObject>
#include <QRegExp>

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
	: VtkVisPipelineItem(algorithm, parentItem,
	                     data), _transformFilter(NULL), _activeAttribute("")
{
	VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
	if (parentItem->parentItem())
	{
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
	: VtkVisPipelineItem(compositeFilter, parentItem,
	                     data), _transformFilter(NULL), _activeAttribute("")
{
}

VtkVisPointSetItem::~VtkVisPointSetItem()
{
	_transformFilter->Delete();
}

void VtkVisPointSetItem::Initialize(vtkRenderer* renderer)
{
	_activeAttribute = "";
	_transformFilter = vtkTransformFilter::New();
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Identity();
	_transformFilter->SetTransform(transform);

	_transformFilter->SetInputConnection(_algorithm->GetOutputPort());
	_transformFilter->Update();

	_renderer = renderer;
	_mapper = QVtkDataSetMapper::New();
	_mapper->InterpolateScalarsBeforeMappingOff();

	// Use a special vtkImageActor instead of vtkActor
	_mapper->SetInputConnection(_transformFilter->GetOutputPort());
	_actor = vtkActor::New();
	static_cast<vtkActor*>(_actor)->SetMapper(_mapper);
	_renderer->AddActor(_actor);

	// Set pre-set properties
	VtkAlgorithmProperties* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(_algorithm);
	if (vtkProps)
		setVtkProperties(vtkProps);

	// Copy properties from parent
	else
	{
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
				parentItem =
				        dynamic_cast<VtkVisPipelineItem*>(parentItem->parentItem());
		}
	}

	// Set active scalar to the desired one from VtkAlgorithmProperties
	// or to match those of the parent.
	if (vtkProps)
	{
		if (vtkProps->GetActiveAttribute().length() > 0)
			this->SetActiveAttribute(vtkProps->GetActiveAttribute());
		else
		{
			VtkVisPointSetItem* visParentItem =
			        dynamic_cast<VtkVisPointSetItem*>(this->parentItem());
			if (visParentItem)
				this->SetActiveAttribute(visParentItem->GetActiveAttribute());
			if (vtkProps->GetTexture() != NULL)
				this->SetActiveAttribute("Solid Color");
		}
	}
}

void VtkVisPointSetItem::setVtkProperties(VtkAlgorithmProperties* vtkProps)
{
	QObject::connect(vtkProps, SIGNAL(ScalarVisibilityChanged(bool)),
	                 _mapper, SLOT(SetScalarVisibility(bool)));

	this->setLookupTableForActiveScalar();

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
	bool onPointData = true;
	if (name.contains(QRegExp("^P-")))
		onPointData = true;
	else if (name.contains(QRegExp("^C-")))
		onPointData = false;
	else if (name.contains("Solid Color"))
	{
		_activeAttribute = "Solid Color";
		_mapper->ScalarVisibilityOff();
		return;
	}
	else
		return;

	// Remove type identifier
	std::string strippedName = QString(name).remove(0, 2).toStdString();
	const char* charName = strippedName.c_str();

	vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->_algorithm->GetOutputDataObject(0));
	if (dataSet)
	{
		if (onPointData)
		{
			vtkPointData* pointData = dataSet->GetPointData();
			if(pointData)
			{
				if(setActiveAttributeOnData(pointData, strippedName))
				{
					_algorithm->SetInputArrayToProcess(
					        0,
					        0,
					        0,
					        vtkDataObject::
					        FIELD_ASSOCIATION_POINTS,
					        charName);
					_mapper->SetScalarModeToUsePointData();
				}
				else
				{
					_activeAttribute = "Solid Color";
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
				if(setActiveAttributeOnData(cellData, strippedName))
				{
					_algorithm->SetInputArrayToProcess(
					        0,
					        0,
					        0,
					        vtkDataObject::
					        FIELD_ASSOCIATION_CELLS,
					        charName);
					_mapper->SetScalarModeToUseCellData();
				}
				else
				{
					_activeAttribute = "Solid Color";
					_mapper->ScalarVisibilityOff();
					return;
				}
			}
		}

		_mapper->SetScalarRange(dataSet->GetScalarRange());
		this->setLookupTableForActiveScalar();
		_mapper->ScalarVisibilityOn();
		//_mapper->Update();	// KR: TODO - this is incredibly slow ... WHY???
		_activeAttribute = name;
	}
}

bool VtkVisPointSetItem::setActiveAttributeOnData(vtkDataSetAttributes* data, std::string& name)
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

void VtkVisPointSetItem::setLookupTableForActiveScalar()
{
	VtkAlgorithmProperties* vtkProps = dynamic_cast<VtkAlgorithmProperties*>(_algorithm);
	if (vtkProps)
	{
		QVtkDataSetMapper* mapper = dynamic_cast<QVtkDataSetMapper*>(_mapper);
		if (mapper)
		{
			if (vtkProps->GetLookupTable(this->GetActiveAttribute()) == NULL) // default color table
			{
				vtkLookupTable* lut = vtkLookupTable::New();
				vtkProps->SetLookUpTable(GetActiveAttribute(), lut);
			}
			else // specific color table

				_mapper->SetLookupTable(vtkProps->GetLookupTable(this->
				                                                 GetActiveAttribute()));

			_mapper->SetScalarRange(_transformFilter->GetOutput()->GetScalarRange());
			//_mapper->Update();  KR: not necessary?!
		}
	}
}

void VtkVisPointSetItem::SetScalarRange(double min, double max)
{
	_mapper->SetScalarRange(min, max);
	_mapper->Update();
}

void VtkVisPointSetItem::setScale(double x, double y, double z) const
{
	if (this->transformFilter())
	{
		vtkTransform* transform =
		        static_cast<vtkTransform*>(this->transformFilter()->GetTransform());
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
		        static_cast<vtkTransform*>(this->transformFilter()->GetTransform());
		double* scale = transform->GetScale();
		transform->Identity();
		transform->Scale(scale);
		transform->Translate(x, y, z);
		this->transformFilter()->Modified();
	}
}

