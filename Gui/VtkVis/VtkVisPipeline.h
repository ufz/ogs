/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-17
 * \brief  Definition of the VtkVisPipeline class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKVISPIPELINE_H
#define VTKVISPIPELINE_H

// ** INCLUDES **
#include "Color.h"
#include "FEMCondition.h"
#include "GeoType.h"
#include "MshEnums.h"
#include "Point.h"
#include "TreeModel.h"

#include <QMap>
#include <QVector>

class vtkAlgorithm;
class vtkDataSet;
class vtkLight;
class vtkPointSet;
class vtkPolyDataAlgorithm;
class vtkRenderer;
class vtkProp3D;
class QModelIndex;
class QString;
class GeoTreeModel;
class ProcessModel;
class MshModel;
class StationTreeModel;
class TreeModel;
class VtkVisPipelineItem;
class VtkMeshSource;

/**
 * \brief VtkVisPipeline manages the VTK visualization.
 * It is a TreeModel and provides functions for adding and removing OGS
 * Model or vtkAlgorithm objects.
 */
class VtkVisPipeline : public TreeModel
{
	Q_OBJECT

public:
	VtkVisPipeline(vtkRenderer* renderer, QObject* parent = 0);

	/// \brief Emits vtkVisPipelineChanged() and calls base class method.
	bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);

	/// \brief Adds a light to the scene at the given coordinates.
	void addLight(const GeoLib::Point &pos);

	/// \brief Returns a light (or NULL) for the given coordinates.
	vtkLight* getLight(const GeoLib::Point &pos) const;

	/// \brief Removes a light at the given coordinates (if possible).
	void removeLight(const GeoLib::Point &pos);

	/// \brief Returns the background-colour of the scene.
	const QColor getBGColor() const;

	/// \brief Sets the background-colour of the scene.
	void setBGColor(const QColor &color);

	/// \brief Returns the QModelIndex of VtkVisPipelineItem which actor
	/// is the given one.
	QModelIndex getIndex(vtkProp3D* actor);

	Qt::ItemFlags flags( const QModelIndex &index ) const;

	/// \brief Loads a vtk object from the given file and adds it to the pipeline.
	void loadFromFile(QString filename);

	/// \brief Defaults to on.
	void resetCameraOnAddOrRemove(bool reset) { _resetCameraOnAddOrRemove = reset; }

	/// \brief Sets a global superelevation factor on all source items and resets
	/// the factor on other items to 1.
	void setGlobalSuperelevation(double factor) const;

	/// \brief Enables / disables backface culling on all actors.
	void setGlobalBackfaceCulling(bool enable) const;

public slots:
	/// \brief Adds the given Model to the pipeline.
	void addPipelineItem(MshModel* model, const QModelIndex &idx);
	void addPipelineItem(GeoTreeModel* model, const std::string &name, GeoLib::GEOTYPE type);
	void addPipelineItem(ProcessModel* model, const FiniteElement::ProcessType pcs_type, FEMCondition::CondType cond_type);
	void addPipelineItem(StationTreeModel* model, const std::string &name);
	QModelIndex addPipelineItem(VtkVisPipelineItem* item, const QModelIndex &parent);

	/// \brief Inserts the vtkAlgorithm as a child of the given QModelIndex to the pipeline.
	QModelIndex addPipelineItem(vtkAlgorithm* source, QModelIndex parent = QModelIndex());

	/// \brief Removes the given Model (and all attached vtkAlgorithms) from the pipeline.
	void removeSourceItem(MshModel* model, const QModelIndex &idx);
	void removeSourceItem(GeoTreeModel* model, const std::string &name, GeoLib::GEOTYPE type);
	void removeSourceItem(ProcessModel* model, const FiniteElement::ProcessType pcs_type, FEMCondition::CondType cond_type);
	void removeSourceItem(StationTreeModel* model, const std::string &name);

	/// \brief Removes the vtkAlgorithm at the given QModelIndex (and all attached
	/// vtkAlgorithms) from the pipeline.
	void removePipelineItem(QModelIndex index);

	/// Checks the quality of a mesh and cal a filter to highlight deformed elements.
	void checkMeshQuality(VtkMeshSource* mesh, MshQualityType::type t);

	void highlightGeoObject(const vtkPolyDataAlgorithm* source, int index);
	void removeHighlightedGeoObject();

private:
	void listArrays(vtkDataSet* dataSet);

	vtkRenderer* _renderer;
	QVector<vtkAlgorithm*> _sources;
	std::list<vtkLight*> _lights;
	QMap<vtkProp3D*, QModelIndex> _actorMap;
	bool _resetCameraOnAddOrRemove;

	QModelIndex _highlighted_geo_index;


signals:
	/// \brief Is emitted when a pipeline item was added or removed.
	void vtkVisPipelineChanged() const;
};

#endif // VTKVISPIPELINE_H
