/**
 * \file VtkVisPipelineItem.cpp
 * 17/2/2010 LB Initial implementation
 *
 * Implementation of VtkVisPipelineItem
 */

// ** INCLUDES **
#include "VtkVisPipelineItem.h"
#include "VtkAlgorithmProperties.h"

#include <vtkAlgorithm.h>
#include <vtkPointSet.h>
#include <vtkDataSetMapper.h>
#include "QVtkDataSetMapper.h"
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkTransformFilter.h>
#include <vtkTransform.h>
#include <vtkTextureMapToPlane.h>

#include <vtkGenericDataObjectWriter.h>

#include <QMessageBox>
#include <QObject>

// export test
#include <vtkPolyDataAlgorithm.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageActor.h>
#include <vtkImageAlgorithm.h>
#include <vtkTubeFilter.h>
#include <vtkTriangleFilter.h>

#include "VtkCompositeFilter.h"

#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>

#ifdef OGS_USE_OPENSG
#include "vtkOsgConverter.h"
#include <OpenSG/OSGSceneFileHandler.h>
#endif

VtkVisPipelineItem::VtkVisPipelineItem(
	vtkAlgorithm* algorithm, TreeItem* parentItem,
	const QList<QVariant> data /*= QList<QVariant>()*/)
: TreeItem(data, parentItem),	_actor(NULL), _algorithm(algorithm), _mapper(NULL), _renderer(NULL),
	  _compositeFilter(NULL), _transformFilter(NULL), _activeAttribute("")
{
	VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
	if (parentItem->parentItem())
	{
		if (dynamic_cast<vtkImageAlgorithm*>(visParentItem->algorithm()))
			_algorithm->SetInputConnection(visParentItem->algorithm()->GetOutputPort());
		else
			_algorithm->SetInputConnection(visParentItem->transformFilter()->GetOutputPort());
	}
}

VtkVisPipelineItem::VtkVisPipelineItem(
	VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
	const QList<QVariant> data /*= QList<QVariant>()*/)
: TreeItem(data, parentItem), 	_actor(NULL), _mapper(NULL), _renderer(NULL),
	  _compositeFilter(compositeFilter), _transformFilter(NULL), _activeAttribute("")
{
	_algorithm = _compositeFilter->GetOutputAlgorithm();
}


VtkVisPipelineItem::~VtkVisPipelineItem()
{
	_renderer->RemoveActor(_actor);
	_mapper->Delete();
	_actor->Delete();
	//_algorithm->Delete();	// TODO: not calling delete causes memoryleak in some cases (e.g. building mesh from images),
	// always calling it causes error when closing program
	delete _compositeFilter;
	if (_transformFilter) _transformFilter->Delete();
}

VtkVisPipelineItem* VtkVisPipelineItem::child( int row ) const
{
	TreeItem* treeItem = TreeItem::child(row);
	if (treeItem)
		return dynamic_cast<VtkVisPipelineItem*>(treeItem);
	else
		return NULL;
}
QVariant VtkVisPipelineItem::data( int column ) const
{
	if (column == 1)
	{
		return isVisible();
	}
	else
		return TreeItem::data(column);
}

bool VtkVisPipelineItem::setData( int column, const QVariant &value )
{
	if (column == 1)
	{
		setVisible(value.toBool());
		return true;
	}
	else
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

void VtkVisPipelineItem::Initialize(vtkRenderer* renderer)
{
	_activeAttribute = "";

	vtkImageAlgorithm* imageAlgorithm = dynamic_cast<vtkImageAlgorithm*>(_algorithm);

	if (imageAlgorithm==NULL) // if algorithm is no image
	{
		_transformFilter = vtkTransformFilter::New();
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		transform->Identity();
		_transformFilter->SetTransform(transform);

		_transformFilter->SetInputConnection(_algorithm->GetOutputPort());
		_transformFilter->Update();
	}

	_renderer = renderer;
	_mapper = QVtkDataSetMapper::New();
	_mapper->InterpolateScalarsBeforeMappingOff();

	// Use a special vtkImageActor instead of vtkActor
	if (imageAlgorithm)
	{
		vtkImageAlgorithm* imageAlg = static_cast<vtkImageAlgorithm*>(_algorithm);
		vtkImageActor* imageActor = vtkImageActor::New();
		imageActor->SetInput(imageAlg->GetOutput());
		_actor = imageActor;
	}
	else
	{
		_mapper->SetInputConnection(_transformFilter->GetOutputPort());
		_actor = vtkActor::New();
		static_cast<vtkActor*>(_actor)->SetMapper(_mapper);
	}
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
	}

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

void VtkVisPipelineItem::setVtkProperties(VtkAlgorithmProperties* vtkProps)
{
	QObject::connect(vtkProps, SIGNAL(ScalarVisibilityChanged(bool)),
		_mapper, SLOT(SetScalarVisibility(bool)));

	vtkImageAlgorithm* imageAlgorithm = dynamic_cast<vtkImageAlgorithm*>(_algorithm);
	if (imageAlgorithm==NULL)
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

int VtkVisPipelineItem::writeToFile(const std::string &filename) const
{
	if (!filename.empty())
	{
		if (filename.substr(filename.size() - 4).find("os") != std::string::npos)
		{
			#ifdef OGS_USE_OPENSG
			vtkOsgConverter osgConverter(static_cast<vtkActor*>(_actor));
			if(osgConverter.WriteAnActor())
				OSG::SceneFileHandler::the().write(osgConverter.GetOsgNode(), filename.c_str());
			#else
			QMessageBox::warning(NULL, "Functionality not implemented",
				"Sorry but this program was not compiled with OpenSG support.");
			#endif
			return 0;
		}

		vtkAlgorithm* alg = this->algorithm();
		vtkPolyDataAlgorithm* algPD = dynamic_cast<vtkPolyDataAlgorithm*>(alg);
		vtkUnstructuredGridAlgorithm* algUG = dynamic_cast<vtkUnstructuredGridAlgorithm*>(alg);
		vtkImageAlgorithm* algI = dynamic_cast<vtkImageAlgorithm*>(alg);
		if (algPD)
		{
//			vtkGenericDataObjectWriter* pdWriter = vtkGenericDataObjectWriter::New();
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
		else if (algI)
		{
			vtkSmartPointer<vtkXMLImageDataWriter> iWriter =
				vtkSmartPointer<vtkXMLImageDataWriter>::New();
			iWriter->SetInput(algI->GetOutputDataObject(0));
			std::string filenameWithExt = filename;
			filenameWithExt.append(".vti");
			iWriter->SetFileName(filenameWithExt.c_str());
			return iWriter->Write();
		}
		std::cout << "VtkVisPipelineItem::writeToFile() - Unknown data type..." << std::endl;
	}
	return 0;
}

vtkProp3D* VtkVisPipelineItem::actor() const
{
	return _actor;
}

void VtkVisPipelineItem::SetScalarVisibility( bool on )
{
	_mapper->SetScalarVisibility(on);
}

void VtkVisPipelineItem::SetActiveAttribute( const QString& name )
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
					_algorithm->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, charName);
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
					_algorithm->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS, charName);
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

bool VtkVisPipelineItem::setActiveAttributeOnData(vtkDataSetAttributes* data, std::string& name)
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

void VtkVisPipelineItem::setLookupTableForActiveScalar()
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
			{
				_mapper->SetLookupTable(vtkProps->GetLookupTable(this->GetActiveAttribute()));
			}

			_mapper->SetScalarRange(_transformFilter->GetOutput()->GetScalarRange());
			//_mapper->Update();  KR: not necessary?!
		}
	}
}

void VtkVisPipelineItem::SetScalarRange(double min, double max)
{
	_mapper->SetScalarRange(min, max);
	_mapper->Update();
}

void VtkVisPipelineItem::setScale(double x, double y, double z) const
{
	if (this->transformFilter())
	{
		vtkTransform* transform =
			static_cast<vtkTransform*>(this->transformFilter()->GetTransform());
		double* trans = transform->GetPosition();
		transform->Identity();
		transform->Scale(x, y, z);
		transform->Translate(trans[0]/x, trans[1]/y, trans[2]/z);
		this->transformFilter()->Modified();
	}

}

void VtkVisPipelineItem::setTranslation(double x, double y, double z) const
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

void VtkVisPipelineItem::setScaleOnChildren(double x, double y, double z) const
{
	for (int i = 0; i < this->childCount(); ++i)
	{
		VtkVisPipelineItem* child = this->child(i);
		child->setScale(x, y, z);
	}
}
