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

#include <QMessageBox>
#include <QObject>

#ifdef OGS_USE_OPENSG
#include <OpenSG/OSGSceneFileHandler.h>
#include <OpenSG/OSGCoredNodePtr.h>
#include <OpenSG/OSGGroup.h>
#include <OpenSG/OSGNode.h>
#include "vtkOsgActor.h"
#endif

// export test
#include <vtkPolyDataAlgorithm.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageActor.h>
#include <vtkImageAlgorithm.h>

#include "VtkCompositeFilter.h"

#include <vtkPointData.h>
#include <vtkCellData.h>

#ifdef OGS_USE_OPENSG
OSG::NodePtr VtkVisPipelineItem::rootNode = NullFC;

	VtkVisPipelineItem::VtkVisPipelineItem(
		vtkAlgorithm* algorithm,
		TreeItem* parentItem,
		const QList<QVariant> data /*= QList<QVariant>()*/)
	: TreeItem(data, parentItem), _algorithm(algorithm), _compositeFilter(NULL)
	{
		VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
		if (visParentItem)
		{
			_parentNode = static_cast<vtkOsgActor*>(visParentItem->actor())->GetOsgRoot();

			if (parentItem->parentItem())
			{
				if (dynamic_cast<vtkImageAlgorithm*>(visParentItem->algorithm()))
					_algorithm->SetInputConnection(visParentItem->algorithm()->GetOutputPort());
				else
					_algorithm->SetInputConnection(visParentItem->transformFilter()->GetOutputPort());
			}
		}
		else
			_parentNode = rootNode;
	}

	VtkVisPipelineItem::VtkVisPipelineItem(
		VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
		const QList<QVariant> data /*= QList<QVariant>()*/ )
		: TreeItem(data, parentItem), _compositeFilter(compositeFilter)
	{
		_algorithm = _compositeFilter->GetOutputAlgorithm();
		VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
		if (visParentItem)
			_parentNode = static_cast<vtkOsgActor*>(visParentItem->actor())->GetOsgRoot();
		else
			_parentNode = rootNode;
	}

#else // OGS_USE_OPENSG
	VtkVisPipelineItem::VtkVisPipelineItem(
		vtkAlgorithm* algorithm, TreeItem* parentItem,
		const QList<QVariant> data /*= QList<QVariant>()*/)
	: TreeItem(data, parentItem), _algorithm(algorithm), _compositeFilter(NULL)
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
	: TreeItem(data, parentItem), _compositeFilter(compositeFilter)
	{
		_algorithm = _compositeFilter->GetOutputAlgorithm();
	}
#endif // OGS_USE_OPENSG


VtkVisPipelineItem::~VtkVisPipelineItem()
{
	#ifdef OGS_USE_OPENSG
		vtkOsgActor* osgActor = dynamic_cast<vtkOsgActor*>(actor());
		if(_parentNode != NullFC)
		{
			OSG::beginEditCP(_parentNode);{
				OSG::RefPtr<OSG::NodePtr> node;
				node = osgActor->GetOsgRoot();
				_parentNode->subChild(node);
			};OSG::endEditCP(_parentNode);
		}
		_renderer->RemoveActor(osgActor);
		_mapper->Delete();
		osgActor->Delete();
	#else // OGS_USE_OPENSG
		_renderer->RemoveActor(_actor);
		_mapper->Delete();
		_actor->Delete();
		//_algorithm->Delete();	// TODO: not calling delete causes memoryleak in some cases (e.g. building mesh from images),
		                        // always calling it causes error when closing program
	#endif // OGS_USE_OPENSG
		delete _compositeFilter;
		_transformFilter->Delete();
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
	_transformFilter = vtkTransformFilter::New();
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Identity();
	_transformFilter->SetTransform(transform);

	_transformFilter->SetInputConnection(_algorithm->GetOutputPort());
	_transformFilter->Update();
	_renderer = renderer;
	_mapper = QVtkDataSetMapper::New();
	_mapper->InterpolateScalarsBeforeMappingOff();


	vtkImageAlgorithm* imageAlgorithm = dynamic_cast<vtkImageAlgorithm*>(_algorithm);

#ifdef OGS_USE_OPENSG
	_actor = vtkOsgActor::New();

	// Transform vtkImageData to vtkPolyData for OpenSG-conversion
	if (imageAlgorithm)
	{
		vtkSmartPointer<vtkTextureMapToPlane> toPolyData =
			vtkSmartPointer<vtkTextureMapToPlane>::New();
		toPolyData->SetInputConnection(imageAlgorithm->GetOutputPort());
		_mapper->SetInputConnection(toPolyData->GetOutputPort());
	}
	else
	{
		_mapper->SetInputConnection(_transformFilter->GetOutputPort());
	}
	_actor->SetMapper(_mapper);

	OSG::beginEditCP(_parentNode);{
		_parentNode->addChild(_actor->GetOsgRoot());
	};OSG::endEditCP(_parentNode);
#else
	_mapper->SetInputConnection(_transformFilter->GetOutputPort());

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
		_actor = vtkActor::New();
		static_cast<vtkActor*>(_actor)->SetMapper(_mapper);
	}
#endif
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
				parentItem = NULL;
			}
			else
				parentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem->parentItem());
		}
	}

}

void VtkVisPipelineItem::setVtkProperties(VtkAlgorithmProperties* vtkProps)
{
	QObject::connect(vtkProps, SIGNAL(ScalarVisibilityChanged(bool)),
		_mapper, SLOT(SetScalarVisibility(bool)));

	//vtkProps->SetLookUpTable("c:/Project/BoreholeColourReferenceMesh.txt"); //HACK ... needs to be put in GUI

	QVtkDataSetMapper* mapper = dynamic_cast<QVtkDataSetMapper*>(_mapper);
	if (mapper)
	{
		if (vtkProps->GetLookupTable() == NULL) // default color table
		{
			vtkLookupTable* lut = vtkLookupTable::New();
			vtkProps->SetLookUpTable(lut);
		}
		else // specific color table
		{
			_mapper->SetLookupTable(vtkProps->GetLookupTable());
		}
		_mapper->SetScalarRange(_transformFilter->GetOutput()->GetScalarRange());
		_mapper->Update();
	}

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
			vtkOsgActor* osgActor = static_cast<vtkOsgActor*>(_actor);
			osgActor->SetVerbose(true);
			osgActor->UpdateOsg();
			OSG::SceneFileHandler::the().write(osgActor->GetOsgRoot(), filename.c_str());
			osgActor->ClearOsg();
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
			vtkSmartPointer<vtkXMLPolyDataWriter> pdWriter =
				vtkSmartPointer<vtkXMLPolyDataWriter>::New();
			pdWriter->SetInput(algPD->GetOutputDataObject(0));
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

void VtkVisPipelineItem::SetActiveAttribute( int arrayIndex, int attributeType )
{
	if (arrayIndex<0)
	{
		return;
	}

	vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->_algorithm->GetOutputDataObject(0));
	bool onPointData(true);
	double* range(NULL);

	int nPointArrays = dataSet->GetPointData()->GetNumberOfArrays();
	int nCellArrays  = dataSet->GetCellData()->GetNumberOfArrays();

	if (arrayIndex==(nPointArrays+nCellArrays))
	{
		_mapper->ScalarVisibilityOff();
		
	}

	//if (nPointArrays+nCellArrays > 0)
	//{
	else 
	{
		if ( arrayIndex > nPointArrays-1 )
		{
			onPointData = false;
			arrayIndex-=nPointArrays;
		}

		if (dataSet)
		{
			if (onPointData)
			{
				vtkPointData* pointData = dataSet->GetPointData();
				const char* charName = pointData->GetArrayName(arrayIndex);
				pointData->SetActiveAttribute(charName, attributeType);
				range = pointData->GetArray(charName)->GetRange();
				_mapper->SetScalarModeToUsePointData();
				_mapper->SetScalarRange(dataSet->GetScalarRange());
			}
			else
			{
				vtkCellData* cellData = dataSet->GetCellData();
				const char* charName = cellData->GetArrayName(arrayIndex);
				cellData->SetActiveAttribute(charName, attributeType);
				range = cellData->GetArray(charName)->GetRange();
				_mapper->SetScalarModeToUseCellData();
				_mapper->SetScalarRange(dataSet->GetScalarRange());
			}

			_mapper->SetScalarRange(range);
			_mapper->ScalarVisibilityOn();
			//_mapper->Update();
		}
	}
		_mapper->Update();
}
