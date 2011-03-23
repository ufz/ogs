/**
 * \file VtkVisPipeline.h
 * 17/2/2010 LB Initial implementation
 *
 */


#ifndef VTKVISPIPELINE_H
#define VTKVISPIPELINE_H

// ** INCLUDES **
#include "Configure.h"
#include "TreeModel.h"
#include "Color.h"
#include "Point.h"
#include "GeoType.h"
#include "MSHEnums.h"

#include <QVector>
#include <QMap>

#ifdef OGS_USE_OPENSG
	#include <OpenSG/OSGSimpleSceneManager.h>
#endif // OGS_USE_OPENSG


class vtkAlgorithm;
class vtkDataSet;
class vtkLight;
class vtkPointSet;
class vtkRenderer;
class vtkProp3D;
class MshModel;
class QModelIndex;
class QString;
class GeoTreeModel;
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

	/// \brief Constructor
#ifdef OGS_USE_OPENSG
	VtkVisPipeline(vtkRenderer* renderer, OSG::SimpleSceneManager* manager, QObject* parent = 0);
#else // OGS_USE_OPENSG
	VtkVisPipeline(vtkRenderer* renderer, QObject* parent = 0);
#endif // OGS_USE_OPENSG

	/// \brief Emits vtkVisPipelineChanged() and calls base class method.
	bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);

	/// \brief Adds a light to the scene at the given coordinates.
	void addLight(const GEOLIB::Point &pos);

	/// \brief Returns a light (or NULL) for the given coordinates.
	vtkLight* getLight(const GEOLIB::Point &pos) const;

	/// \brief Removes a light at the given coordinates (if possible).
	void removeLight(const GEOLIB::Point &pos);

	/// \brief Returns the background-colour of the scene.
	const QColor getBGColor() const;

	/// \brief Sets the background-colour of the scene.
	void setBGColor(const QColor &color);

	/// \brief Returns the QModelIndex of VtkVisPipelineItem which actor
	/// is the given one.
	QModelIndex getIndex(vtkProp3D* actor);

	Qt::ItemFlags flags( const QModelIndex &index ) const;

	/// @brief Loads a vtk object from the given file and adds it to the pipeline.
	void loadFromFile(QString filename);

	///
	void resetCameraOnAddOrRemove(bool reset) { _resetCameraOnAddOrRemove = reset; }

public slots:
	/// \brief Adds the given Model to the pipeline.
	void addPipelineItem(MshModel* model, const QModelIndex &idx);
	void addPipelineItem(GeoTreeModel* model, const std::string &name, GEOLIB::GEOTYPE type);
	void addPipelineItem(StationTreeModel* model, const std::string &name);
	void addPipelineItem(VtkVisPipelineItem* item, const QModelIndex &parent);

	/// \brief Inserts the vtkAlgorithm as a child of the given QModelIndex to the pipeline.
	void addPipelineItem(vtkAlgorithm* source, QModelIndex parent = QModelIndex());

	/// \brief Removes the given Model (and all attached vtkAlgorithms) from the pipeline.
	void removeSourceItem(MshModel* model, const QModelIndex &idx);
	void removeSourceItem(GeoTreeModel* model, const std::string &name, GEOLIB::GEOTYPE type);
	void removeSourceItem(StationTreeModel* model, const std::string &name);

	/// \brief Removes the vtkAlgorithm at the given QModelIndex (and all attached
	/// vtkAlgorithms) from the pipeline.
	void removePipelineItem(QModelIndex index);

	/// Checks the quality of a mesh and cal a filter to highlight deformed elements.
	void checkMeshQuality(VtkMeshSource* mesh, MshQualityType::type t);

private:
	void listArrays(vtkDataSet* dataSet);

	vtkRenderer* _renderer;
	QVector<vtkAlgorithm*> _sources;
	std::list<vtkLight*> _lights;
	QMap<vtkProp3D*, QModelIndex> _actorMap;
	bool _resetCameraOnAddOrRemove;

#ifdef OGS_USE_OPENSG
	OSG::SimpleSceneManager* _sceneManager;
#endif // OGS_USE_OPENSG



signals:
	/// \brief Is emitted when a pipeline item was added or removed.
	void vtkVisPipelineChanged();

};

#endif // VTKVISPIPELINE_H
