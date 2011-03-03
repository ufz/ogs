/**
 * \file VtkVisPipeline.cpp
 * 17/2/2010 LB Initial implementation
 *
 * Implementation of VtkVisPipeline
 */

// ** INCLUDES **
#include "VtkVisPipeline.h"

//#include "Model.h"
#include "TreeModel.h"
#include "MshModel.h"
#include "MshItem.h"
#include "GeoTreeModel.h"
#include "StationTreeModel.h"
#include "VtkVisPipelineItem.h"
#include "VtkMeshSource.h"
#include "VtkAlgorithmProperties.h"
#include "VtkTrackedCamera.h"

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkAlgorithm.h>
#include <vtkPointSet.h>
#include <vtkProp3D.h>
#include <vtkLight.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkImageReader2.h>
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <QString>
#include <QTime>
#include <QFileInfo>
#include <QColor>
#include <QSettings>

#ifdef OGS_USE_OPENSG
#include "vtkOsgActor.h"
VtkVisPipeline::VtkVisPipeline(vtkRenderer* renderer, OSG::SimpleSceneManager* manager, QObject* parent /*= 0*/)
: TreeModel(parent), _renderer(renderer), _sceneManager(manager)
{
	QList<QVariant> rootData;
	rootData << "Object name" << "Visible";
	delete _rootItem;
	_rootItem = new TreeItem(rootData, NULL);
	VtkVisPipelineItem::rootNode = _sceneManager->getRoot();

	QSettings settings("UFZ", "OpenGeoSys-5");
	QVariant backgroundColorVariant = settings.value("VtkBackgroundColor");
	if (backgroundColorVariant != QVariant())
		this->setBGColor(backgroundColorVariant.value<QColor>());
}
#else // OGS_USE_OPENSG
VtkVisPipeline::VtkVisPipeline( vtkRenderer* renderer, QObject* parent /*= 0*/ )
: TreeModel(parent), _renderer(renderer)
{
	QList<QVariant> rootData;
	rootData << "Object name" << "Visible";
	delete _rootItem;
	_rootItem = new TreeItem(rootData, NULL);
	//_renderer->SetBackground(1,1,1);
}
#endif // OGS_USE_OPENSG

bool VtkVisPipeline::setData( const QModelIndex &index, const QVariant &value,
	int role /* = Qt::EditRole */ )
{
	emit vtkVisPipelineChanged();

	return TreeModel::setData(index, value, role);
}

void VtkVisPipeline::addLight(const GEOLIB::Point &pos)
{
	double lightPos[3];
	for (std::list<vtkLight*>::iterator it = _lights.begin(); it != _lights.end(); ++it)
	{
		(*it)->GetPosition(lightPos);
		if (pos[0] == lightPos[0] && pos[1] == lightPos[1] && pos[2] == lightPos[2]) return;
	}
	vtkLight* l = vtkLight::New();
	l->SetPosition(pos[0], pos[1], pos[2]);
	_renderer->AddLight(l);
	_lights.push_back(l);
}

vtkLight* VtkVisPipeline::getLight(const GEOLIB::Point &pos) const
{
	double lightPos[3];
	for (std::list<vtkLight*>::const_iterator it = _lights.begin(); it != _lights.end(); ++it)
	{
		(*it)->GetPosition(lightPos);
		if (pos[0] == lightPos[0] && pos[1] == lightPos[1] && pos[2] == lightPos[2]) return (*it);
	}
	return NULL;
}

void VtkVisPipeline::removeLight(const GEOLIB::Point &pos)
{
	double lightPos[3];
	for (std::list<vtkLight*>::iterator it = _lights.begin(); it != _lights.end(); ++it)
	{
		(*it)->GetPosition(lightPos);
		if (pos[0] == lightPos[0] && pos[1] == lightPos[1] && pos[2] == lightPos[2])
		{
			_renderer->RemoveLight(*it);
			(*it)->Delete();
			_lights.erase(it);
			return;
		}
	}
}

const QColor VtkVisPipeline::getBGColor() const
{
	double* color = _renderer->GetBackground();
	QColor c(static_cast<int>(color[0]*255), static_cast<int>(color[1]*255), static_cast<int>(color[2]*255));
	return c;
}

void VtkVisPipeline::setBGColor(const QColor &color)
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	settings.setValue("VtkBackgroundColor", color);
	_renderer->SetBackground(color.redF(), color.greenF(), color.blueF());
}

QModelIndex VtkVisPipeline::getIndex( vtkProp3D* actor )
{
	return _actorMap.value(actor, QModelIndex());
}

Qt::ItemFlags VtkVisPipeline::flags( const QModelIndex &index ) const
{
	Qt::ItemFlags defaultFlags = Qt::ItemIsEnabled | Qt::ItemIsSelectable;

	if (!index.isValid())
		return Qt::ItemIsEnabled;

	//if (index.column() == 1)
	//	defaultFlags |= Qt::ItemIsEditable;

	return defaultFlags;
}

void VtkVisPipeline::loadFromFile(QString filename)
{
	#ifndef NDEBUG
	    	 QTime myTimer;
	    	 myTimer.start();
			std::cout << "VTK Read: " << filename.toStdString() <<
				std::endl << std::flush;
	#endif

	if (filename.size() > 0)
	{
		vtkXMLDataReader* reader;
		if (filename.endsWith("vti"))
			reader = vtkXMLImageDataReader::New();
		else if (filename.endsWith("vtr"))
			reader = vtkXMLRectilinearGridReader::New();
		else if (filename.endsWith("vts"))
			reader = vtkXMLStructuredGridReader::New();
		else if (filename.endsWith("vtp"))
			reader = vtkXMLPolyDataReader::New();
		else if (filename.endsWith("vtu"))
			reader = vtkXMLUnstructuredGridReader::New();
		else if (filename.endsWith("vtk"))
		{
			vtkGenericDataObjectReader* oldStyleReader = vtkGenericDataObjectReader::New();
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
				std::cout << "Error loading vtk file: not a valid vtkDataSet." << std::endl;

			return;
		}
		else
			return;

		reader->SetFileName(filename.toStdString().c_str());
		// TODO: insert ReadAllScalarsOn()-equivalent for xml-file-reader here, otherwise arrays are not available in GUI!
		reader->Update();
		//std::cout << "#cell scalars: " << reader->GetNumberOfCellArrays() << std::endl;
		//std::cout << "#point scalars: " << reader->GetNumberOfPointArrays() << std::endl;

		vtkDataSet* dataSet = reader->GetOutputAsDataSet();
		if (dataSet)
		{
			this->listArrays(dataSet);
			addPipelineItem(reader);
		}
		else
			std::cout << "Error loading vtk file: not a valid vtkDataSet." << std::endl;

		//reader->Delete();
	}

	#ifndef NDEBUG
	    	 std::cout << myTimer.elapsed() << " ms" << std::endl;
	#endif
}

void VtkVisPipeline::addPipelineItem(GeoTreeModel* model, const std::string &name, GEOLIB::GEOTYPE type)
{
	addPipelineItem(model->vtkSource(name, type));
}

void VtkVisPipeline::addPipelineItem(StationTreeModel* model, const std::string &name)
{
	addPipelineItem(model->vtkSource(name));
}

void VtkVisPipeline::addPipelineItem(MshModel* model, const QModelIndex &idx)
{
	addPipelineItem(static_cast<MshItem*>(model->getItem(idx))->vtkSource());
}

void VtkVisPipeline::addPipelineItem(VtkVisPipelineItem* item, const QModelIndex &parent)
{
	item->Initialize(_renderer);
	TreeItem* parentItem = item->parentItem();
	parentItem->appendChild(item);

	if (parent.isValid())  // KR scale children according to parent
	{
		double* scale = static_cast<VtkVisPipelineItem*>(parentItem)->actor()->GetScale();
		item->actor()->SetScale(scale);
	}

	int parentChildCount = parentItem->childCount();
	QModelIndex newIndex = index(parentChildCount - 1, 0, parent);

	_renderer->ResetCamera(_renderer->ComputeVisiblePropBounds());
	_actorMap.insert(item->actor(), newIndex);

	// Do not interpolate images
#ifndef OGS_USE_OPENSG
	if (dynamic_cast<vtkImageAlgorithm*>(item->algorithm()))
		static_cast<vtkImageActor*>(item->actor())->InterpolateOff();
#endif // OGS_USE_OPENSG

	reset();
	emit vtkVisPipelineChanged();
}

void VtkVisPipeline::addPipelineItem( vtkAlgorithm* source,
									  QModelIndex parent /* = QModelindex() */)
{
	TreeItem* parentItem = getItem(parent);

	// If the parent is not the root TreeItem
	//if (parent.isValid())
	//	VtkVisPipelineItem* visParentItem = static_cast<VtkVisPipelineItem*>(parentItem);

	QList<QVariant> itemData;
	QString itemName;
	if (!parent.isValid())	// if source object
	{
		vtkGenericDataObjectReader* reader = dynamic_cast<vtkGenericDataObjectReader*>(source);
		vtkImageReader2* imageReader = dynamic_cast<vtkImageReader2*>(source);
		VtkAlgorithmProperties* props = dynamic_cast<VtkAlgorithmProperties*>(source);
		if (reader)
		{
			QFileInfo fi(QString(reader->GetFileName()));
			itemName = fi.fileName();
		}
		else if (imageReader)
		{
			QFileInfo fi(QString(imageReader->GetFileName()));
			itemName = fi.fileName();
		}
		else if (props)
		{
			QFileInfo fi(props->GetName());
			itemName = fi.fileName();
		}
		else
			itemName = QString(source->GetClassName());
	}
	else
		itemName = QString(source->GetClassName());
	itemData << itemName << true;


	VtkVisPipelineItem* item = new VtkVisPipelineItem(source, parentItem, itemData);
	this->addPipelineItem(item, parent);

#ifdef OGS_USE_OPENSG
	_sceneManager->showAll();
#endif // OGS_USE_OPENSG
}

void VtkVisPipeline::removeSourceItem(GeoTreeModel* model, const std::string &name, GEOLIB::GEOTYPE type)
{
	for (int i = 0; i < _rootItem->childCount(); i++)
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
	for (int i = 0; i < _rootItem->childCount(); i++)
	{
		VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(getItem(index(i, 0)));
		if (item->algorithm() == model->vtkSource(name))
		{
			removePipelineItem(index(i, 0));
			return;
		}
	}
}

void VtkVisPipeline::removeSourceItem(MshModel* model, const QModelIndex &idx)
{
	MshItem* sItem = static_cast<MshItem*>(model->getItem(idx));

	for (int i = 0; i < _rootItem->childCount(); i++)
	{
		VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(getItem(index(i, 0)));
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
		return;

	QMap<vtkProp3D*, QModelIndex>::iterator it = _actorMap.begin();
	while (it != _actorMap.end())
	{
		QModelIndex itIndex = it.value();
		if (itIndex == index)
		{
			_actorMap.erase(it);
			break;
		}
		++it;
	}

	//TreeItem* item = getItem(index);
	removeRows(index.row(), 1, index.parent());

	_renderer->ResetCamera(_renderer->ComputeVisiblePropBounds());
	emit vtkVisPipelineChanged();
}

void VtkVisPipeline::listArrays(vtkDataSet* dataSet)
{
	if (dataSet)
	{
		vtkPointData* pointData = dataSet->GetPointData();
		std::cout << "  #point data arrays: " << pointData->GetNumberOfArrays() << std::endl;
		for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
			std::cout << "    Name: " << pointData->GetArrayName(i) << std::endl;

		vtkCellData* cellData = dataSet->GetCellData();
		std::cout << "  #cell data arrays: " << cellData->GetNumberOfArrays() << std::endl;
		for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
			std::cout << "    Name: " << cellData->GetArrayName(i) << std::endl;
	}
	else
		std::cout << "Error loading vtk file: not a valid vtkDataSet." << std::endl;
}