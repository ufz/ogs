/**
 * \file GEOModels.h
 * 9/2/2010 LB Initial implementation
 *
 */


#ifndef GEOMODELS_H
#define GEOMODELS_H

// ** INCLUDES **
#include <QObject>
#include "GEOObjects.h"
#

class Model;
class PntsModel;
class StationTreeModel;
class PolylinesModel;
class SurfaceModel;
class TreeModel;

/**
 * GEOModels glues together the data management class GEOObjects and the GUI.
 * It inherits from GEOLIB::GEOObjects and additionally emits signals when
 * data objects are modified. The GUI connects to these signals. Model instances
 * are created for every data object.
 */
class GEOModels : public QObject, public GEOLIB::GEOObjects
{
	Q_OBJECT

public:
	GEOModels(QObject* parent = 0);
	~GEOModels();

	StationTreeModel* getStationModel() { return _stationModel; }

public slots:
	void addPointVec(std::vector<GEOLIB::Point*> *points, std::string &name, std::map<std::string, size_t>* name_pnt_id_map = NULL);
	bool appendPointVec(const std::vector<GEOLIB::Point*> &points, std::string &name);
	bool removePointVec(const std::string &name);

	void addStationVec(std::vector<GEOLIB::Point*> *stations, std::string &name, const GEOLIB::Color* const color);
	void filterStationVec(const std::string &name, const std::vector<PropertyBounds> &bounds);
	bool removeStationVec(const std::string &name);

	void addPolylineVec(std::vector<GEOLIB::Polyline*> *lines, const std::string &name, std::map<std::string,size_t>* ply_names = NULL);
	bool appendPolylineVec(const std::vector<GEOLIB::Polyline*> &polylines, std::string &name);
	bool removePolylineVec(const std::string &name);

	void addSurfaceVec(std::vector<GEOLIB::Surface*> *surfaces, const std::string &name, std::map<std::string,size_t>* sfc_names = NULL);
	bool removeSurfaceVec(const std::string &name);

protected:
	std::vector<PntsModel*> _pntModels;
	std::vector<PolylinesModel*> _lineModels;
	std::vector<SurfaceModel*> _surfaceModels;;
	StationTreeModel* _stationModel;

private:

signals:
	void pointModelAdded(Model* model);
	void pointModelRemoved(Model* model);

	void stationVectorAdded(StationTreeModel* model, std::string name);
	void stationVectorRemoved(StationTreeModel* model, std::string name);

	void polylineModelAdded(Model* model);
	void polylineModelRemoved(Model* model);

	void surfaceModelAdded(Model* model);
	void surfaceModelRemoved(Model* model);

};

#endif // GEOMODELS_H
