/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-17
 * \brief  Implementation of the VtkVisPipeline class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisPipeline.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// MathLib
#include "InterpolationAlgorithms/LinearIntervalInterpolation.h"

//#include "Model.h"
#include "GeoTreeModel.h"
#include "MeshQuality/AngleSkewMetric.h"
#include "MeshQuality/AreaMetric.h"
#include "MeshQuality/VolumeMetric.h"
#include "MeshQuality/EdgeRatioMetric.h"
#include "MeshQuality/RadiusEdgeRatioMetric.h"
#include "MshItem.h"
#include "MshModel.h"
#include "StationTreeModel.h"
#include "TreeModel.h"
#include "VtkAlgorithmProperties.h"
#include "VtkCompositeNodeSelectionFilter.h"
#include "VtkCompositeElementSelectionFilter.h"
#include "VtkCompositeGeoObjectFilter.h"
#include "VtkFilterFactory.h"
#include "VtkVisImageItem.h"
#include "VtkVisPipelineItem.h"
#include "VtkVisPointSetItem.h"

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

VtkVisPipeline::VtkVisPipeline( vtkRenderer* renderer, QObject* parent /*= 0*/ )
	: TreeModel(parent), _renderer(renderer), _highlighted_geo_index(QModelIndex()), _highlighted_mesh_component(QModelIndex())
{
	QList<QVariant> rootData;
	rootData << "Object name" << "Visible";
	delete _rootItem;
	_rootItem = new TreeItem(rootData, NULL);

	QSettings settings;
	QVariant backgroundColorVariant = settings.value("VtkBackgroundColor");
	if (backgroundColorVariant != QVariant())
		this->setBGColor(backgroundColorVariant.value<QColor>());

	_resetCameraOnAddOrRemove = settings.value("resetViewOnLoad", true).toBool();
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
	for (std::list<vtkLight*>::iterator it = _lights.begin(); it != _lights.end(); ++it)
	{
		(*it)->GetPosition(lightPos);
		if (pos[0] == lightPos[0] && pos[1] == lightPos[1] && pos[2] == lightPos[2])
			return;
	}
	vtkLight* l = vtkLight::New();
	l->SetPosition(pos[0], pos[1], pos[2]);
	_renderer->AddLight(l);
	_lights.push_back(l);
}

vtkLight* VtkVisPipeline::getLight(const GeoLib::Point &pos) const
{
	double lightPos[3];
	for (std::list<vtkLight*>::const_iterator it = _lights.begin(); it != _lights.end(); ++it)
	{
		(*it)->GetPosition(lightPos);
		if (pos[0] == lightPos[0] && pos[1] == lightPos[1] && pos[2] == lightPos[2])
			return *it;
	}
	return NULL;
}

void VtkVisPipeline::removeLight(const GeoLib::Point &pos)
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
	QColor c(static_cast<int>(color[0] * 255),
	         static_cast<int>(color[1] * 255),
	         static_cast<int>(color[2] * 255));
	return c;
}

void VtkVisPipeline::setBGColor(const QColor &color)
{
	QSettings settings;
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
	INFO("VTK Read: %s.", filename.toStdString().c_str());
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
			ERR("VtkVisPipeline::loadFromFile(): not a valid vtkDataSet.");

		//reader->Delete();
	}

#ifndef NDEBUG
	INFO("%d ms", myTimer.elapsed());
#endif
}

void VtkVisPipeline::setGlobalSuperelevation(double factor) const
{
	// iterate over all source items
	for (int i = 0; i < _rootItem->childCount(); ++i)
	{
		VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(_rootItem->child(i));
		item->setScale(1.0, 1.0, factor);

		// recursively set on all child items
		item->setScaleOnChildren(1.0, 1.0, 1.0);
	}

	emit vtkVisPipelineChanged();
}

void VtkVisPipeline::setGlobalBackfaceCulling(bool enable) const
{
	// iterate over all source items
	for (int i = 0; i < _rootItem->childCount(); ++i)
	{
		VtkVisPipelineItem* item = static_cast<VtkVisPipelineItem*>(_rootItem->child(i));
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

void VtkVisPipeline::addPipelineItem(MshModel* model, const QModelIndex &idx)
{
	addPipelineItem(static_cast<MshItem*>(model->getItem(idx))->vtkSource());
}

QModelIndex VtkVisPipeline::addPipelineItem(VtkVisPipelineItem* item, const QModelIndex &parent)
{
	item->Initialize(_renderer);
	TreeItem* parentItem = item->parentItem();
	parentItem->appendChild(item);

	if (!parent.isValid()) // Set global superelevation on source objects
	{
		QSettings settings;
		if (dynamic_cast<vtkImageAlgorithm*>(item->algorithm()) == NULL) // if not an image
			item->setScale(1.0, 1.0, settings.value("globalSuperelevation", 1.0).toDouble());
	}

	int parentChildCount = parentItem->childCount();
	QModelIndex newIndex = index(parentChildCount - 1, 0, parent);

	if (_resetCameraOnAddOrRemove || _rootItem->childCount() == 1)
		_renderer->ResetCamera(_renderer->ComputeVisiblePropBounds());
	_actorMap.insert(item->actor(), newIndex);

	// Do not interpolate images
	if (dynamic_cast<vtkImageAlgorithm*>(item->algorithm()))
		static_cast<vtkImageActor*>(item->actor())->InterpolateOff();

	reset();
	emit vtkVisPipelineChanged();
	emit itemSelected(newIndex);

	return newIndex;
}

QModelIndex VtkVisPipeline::addPipelineItem( vtkAlgorithm* source, QModelIndex parent /* = QModelindex() */)
{
	QString itemName("");

	if (!parent.isValid()) // if source object
	{
		QFileInfo fi;
		vtkGenericDataObjectReader* old_reader = dynamic_cast<vtkGenericDataObjectReader*>(source);
		vtkXMLReader* new_reader = dynamic_cast<vtkXMLReader*>(source);
		vtkImageReader2* image_reader = dynamic_cast<vtkImageReader2*>(source);
		VtkAlgorithmProperties* props = dynamic_cast<VtkAlgorithmProperties*>(source);
		if (old_reader)
			fi.setFile(QString(old_reader->GetFileName()));
		else if (new_reader)
			fi.setFile(QString(new_reader->GetFileName()));
		else if (image_reader)
			fi.setFile(QString(image_reader->GetFileName()));
		else if (props)
			fi.setFile(props->GetName());

		itemName = fi.fileName();
	}

	if (itemName.isEmpty())	// not sure if this may be true at any time
		itemName = QString(source->GetClassName());

	QList<QVariant> itemData;
	itemData << itemName << true;

	VtkVisPipelineItem* item(NULL);
	if (dynamic_cast<vtkImageAlgorithm*>(source))
		item = new VtkVisImageItem(source, getItem(parent), itemData);
	else
		item = new VtkVisPointSetItem(source, getItem(parent), itemData);
	return this->addPipelineItem(item, parent);
}

void VtkVisPipeline::removeSourceItem(GeoTreeModel* model,
                                      const std::string &name,
                                      GeoLib::GEOTYPE type)
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

	if (_resetCameraOnAddOrRemove)
		_renderer->ResetCamera(_renderer->ComputeVisiblePropBounds());
	emit vtkVisPipelineChanged();
}

void VtkVisPipeline::listArrays(vtkDataSet* dataSet)
{
	if (dataSet)
	{
		vtkPointData* pointData = dataSet->GetPointData();
		INFO("  #point data arrays: %d", pointData->GetNumberOfArrays());
		for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
			INFO("    Name: %s", pointData->GetArrayName(i));

		vtkCellData* cellData = dataSet->GetCellData();
		INFO("  #cell data arrays: %d", cellData->GetNumberOfArrays());
		for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
			INFO("    Name: %s", cellData->GetArrayName(i));
	}
	else
		ERR("VtkVisPipeline::listArrays(): not a valid vtkDataSet.");
}

void VtkVisPipeline::checkMeshQuality(InSituLib::VtkMappedMeshSource* source, MeshQualityType t)
{
	if (source)
	{
		const MeshLib::Mesh* mesh = source->GetMesh();
		MeshLib::ElementQualityMetric* quality_tester (nullptr);
		if (t == MeshQualityType::EDGERATIO)
			quality_tester = new MeshLib::EdgeRatioMetric(*mesh);
		else if (t == MeshQualityType::AREA)
			quality_tester = new MeshLib::AreaMetric(*mesh);
		else if (t == MeshQualityType::VOLUME)
			quality_tester = new MeshLib::VolumeMetric(*mesh);
		else if (t == MeshQualityType::EQUIANGLESKEW)
			quality_tester = new MeshLib::AngleSkewMetric(*mesh);
		else if (t == MeshQualityType::RADIUSEDGERATIO)
			quality_tester = new MeshLib::RadiusEdgeRatioMetric(*mesh);
		else
		{
			ERR("VtkVisPipeline::checkMeshQuality(): Unknown MeshQualityType.");
			delete quality_tester;
			return;
		}
		quality_tester->calculateQuality();

		std::vector<double> quality (quality_tester->getElementQuality());
		// transform area and volume criterion values to [0, 1]
		if (t == MeshQualityType::AREA || t == MeshQualityType::VOLUME) {
			try {
				MathLib::LinearIntervalInterpolation<double> lin_intpol(quality_tester->getMinValue(), quality_tester->getMaxValue(), 0, 1);
				const size_t n_quality(quality.size());
				for (size_t k(0); k<n_quality; k++)
					quality[k] = lin_intpol(quality[k]);
			} catch (std::runtime_error& exception) {
				ERR("VtkVisPipeline::checkMeshQuality(): %s", exception.what());
			}
		}

		int nSources = this->_rootItem->childCount();
		for (int i = 0; i < nSources; i++)
		{
			VtkVisPipelineItem* parentItem =
			        static_cast<VtkVisPipelineItem*>(_rootItem->child(i));
			if (parentItem->algorithm() == source)
			{
				QList<QVariant> itemData;
				itemData << "MeshQuality: " + QString::fromStdString(
				        MeshQualityType2String(t)) << true;

				VtkCompositeFilter* filter =
				        VtkFilterFactory::CreateCompositeFilter(
				                "VtkCompositeElementSelectionFilter",
				                parentItem->transformFilter());
				static_cast<VtkCompositeElementSelectionFilter*>(filter)->setSelectionArray("Selection", quality);
				VtkVisPointSetItem* item = new VtkVisPointSetItem(filter, parentItem, itemData);
				this->addPipelineItem(item, this->createIndex(i, 0, item));
				break;
			}
		}

		// *** construct and write histogram
		// simple suggestion: number of classes with Sturges criterion
		size_t nclasses (static_cast<size_t>(1 + 3.3 * log (static_cast<float>(mesh->getNElements()))));
//			bool ok;
//			size_t size (static_cast<size_t>(QInputDialog::getInt(NULL, "OGS-Histogram", "number of histogram classes/spins (min: 1, max: 10000)", static_cast<int>(nclasses), 1, 10000, 1, &ok)));
//			if (ok) ...

		BaseLib::Histogram<double> histogram (quality_tester->getHistogram(nclasses));
		std::ofstream out ("mesh_histogram.txt");
		if (out) {
			out << "# histogram depicts mesh quality criterion " << MeshQualityType2String(t)
				<< " for mesh " << source->GetMesh()->getName() << "\n";
			nclasses = histogram.getNrBins();
			std::vector<size_t> const& bin_cnts(histogram.getBinCounts());
			const double min (histogram.getMinimum());
			const double bin_width (histogram.getBinWidth());

			for (size_t k(0); k < nclasses; k++)
				out << min+k*bin_width << " " << bin_cnts[k] << "\n";
			out.close ();
		} else {
			std::cerr << "could not open file mesh_histgram.txt\n";
		}

//		size_t size (100);
//		std::vector<size_t> histogram (size,0);
//		checker->getHistogram(histogram);
//		std::ofstream out ("mesh_histogram.txt");
//		const size_t histogram_size (histogram.size());
//		for (size_t k(0); k < histogram_size; k++)
//			out << k / static_cast<double>(histogram_size) << " " << histogram[k] << "\n";
//		out.close ();

		delete quality_tester;
	}
}

void VtkVisPipeline::highlightGeoObject(const vtkPolyDataAlgorithm* source, int index)
{
	this->removeHighlightedGeoObject();
	int nSources = this->_rootItem->childCount();
	for (int i = 0; i < nSources; i++)
	{
		VtkVisPipelineItem* parentItem = static_cast<VtkVisPipelineItem*>(_rootItem->child(i));
		if (parentItem->algorithm() == source)
		{
			QList<QVariant> itemData;
			itemData << "Selected GeoObject" << true;

			VtkCompositeFilter* filter = VtkFilterFactory::CreateCompositeFilter(
															"VtkCompositeGeoObjectFilter",
															parentItem->transformFilter());
			static_cast<VtkCompositeGeoObjectFilter*>(filter)->SetIndex(index);
			VtkVisPointSetItem* item = new VtkVisPointSetItem(filter, parentItem, itemData);
			QModelIndex parent_index = static_cast<TreeModel*>(this)->index(i, 0, QModelIndex());
			_highlighted_geo_index = this->addPipelineItem(item, parent_index);
			break;
		}
	}
}

void VtkVisPipeline::removeHighlightedGeoObject()
{
	if (_highlighted_geo_index != QModelIndex())
	{
		this->removePipelineItem(_highlighted_geo_index);
		_highlighted_geo_index = QModelIndex();
	}
}

void VtkVisPipeline::highlightMeshComponent(vtkUnstructuredGridAlgorithm const*const source, unsigned index, bool is_element)
{
	int nSources = this->_rootItem->childCount();
	for (int i = 0; i < nSources; i++)
	{
		VtkVisPipelineItem* parentItem = static_cast<VtkVisPipelineItem*>(_rootItem->child(i));
		if (parentItem->algorithm() == source)
		{
			QList<QVariant> itemData;
			itemData << "Selected Mesh Component" << true;
			QList<QVariant> selected_index;
			selected_index << index << index;

			VtkCompositeFilter* filter (nullptr);
			if (is_element)
			{
				filter = VtkFilterFactory::CreateCompositeFilter("VtkCompositeElementSelectionFilter",
																  parentItem->transformFilter());
				static_cast<VtkCompositeElementSelectionFilter*>(filter)->setSelectionArray("vtkIdFilter_Ids");
				static_cast<VtkCompositeElementSelectionFilter*>(filter)->SetUserVectorProperty("Threshold Between", selected_index);
			}
			else
			{
				filter = VtkFilterFactory::CreateCompositeFilter("VtkCompositeNodeSelectionFilter",
																  parentItem->transformFilter());
				std::vector<unsigned> indeces(1);
				indeces[0] = index;
				static_cast<VtkCompositeNodeSelectionFilter*>(filter)->setSelectionArray(indeces);
			}
			VtkVisPointSetItem* item = new VtkVisPointSetItem(filter, parentItem, itemData);
			QModelIndex parent_index = static_cast<TreeModel*>(this)->index(i, 0, QModelIndex());
			_highlighted_mesh_component = this->addPipelineItem(item, parent_index);
			break;
		}
	}
}

void VtkVisPipeline::removeHighlightedMeshComponent()
{
	if (_highlighted_mesh_component != QModelIndex())
	{
		this->removePipelineItem(_highlighted_mesh_component);
		_highlighted_mesh_component = QModelIndex();
	}
}
