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
#include "GeoType.h"

class GeoTreeModel;
class Model;
class PntsModel;
class StationTreeModel;
class PolylinesModel;
class SurfaceModel;
class TreeModel;

/**
 * \brief GEOModels connects the data management class GEOObjects and the GUI.
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

	GeoTreeModel* getGeoModel() { return _geoModel; }
	StationTreeModel* getStationModel() { return _stationModel; }

public slots:
	/// Removes all parts (points, lines, surfaces) of the geometry with the given name.
	void removeGeometry(std::string geo_name, GEOLIB::GEOTYPE type);

	void addPointVec(std::vector<GEOLIB::Point*> *points, std::string &name, std::map<std::string, size_t>* name_pnt_id_map = NULL);
	bool appendPointVec(const std::vector<GEOLIB::Point*> &points, std::string &name);
	bool removePointVec(const std::string &name);

	void addStationVec(std::vector<GEOLIB::Point*> *stations, std::string &name, const GEOLIB::Color* const color);
	void filterStationVec(const std::string &name, const std::vector<PropertyBounds> &bounds);
	bool removeStationVec(const std::string &name);

	void addPolylineVec(std::vector<GEOLIB::Polyline*> *lines, const std::string &name, std::map<std::string,size_t>* ply_names = NULL);
	bool appendPolylineVec(const std::vector<GEOLIB::Polyline*> &polylines, const std::string &name);
	bool removePolylineVec(const std::string &name);

	void addSurfaceVec(std::vector<GEOLIB::Surface*> *surfaces, const std::string &name, std::map<std::string,size_t>* sfc_names = NULL);
	bool appendSurfaceVec(const std::vector<GEOLIB::Surface*> &surfaces, const std::string &name);
	bool removeSurfaceVec(const std::string &name);

	/// Calls all necessary functions to connect polyline-segments and update all views and windows.
	void connectPolylineSegments(const std::string &geoName, std::vector<size_t> indexlist, bool closePly, bool triangulatePly);


protected:
	GeoTreeModel* _geoModel;
	StationTreeModel* _stationModel;

private:

signals:
	void geoDataAdded(GeoTreeModel*, std::string, GEOLIB::GEOTYPE);
	void geoDataRemoved(GeoTreeModel*, std::string, GEOLIB::GEOTYPE);

	void stationVectorAdded(StationTreeModel* model, std::string name);
	void stationVectorRemoved(StationTreeModel* model, std::string name);
};

#endif // GEOMODELS_H
