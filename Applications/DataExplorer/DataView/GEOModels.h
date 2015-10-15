/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-09
 * \brief  Definition of the GEOModels class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOMODELS_H
#define GEOMODELS_H

#include <QObject>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/GeoType.h"

class GeoTreeModel;
class Model;
class PntsModel;
class StationTreeModel;
class PolylinesModel;
class SurfaceModel;
class TreeModel;

/**
 * \brief GEOModels connects the data management class GEOObjects and the GUI.
 * It inherits from GeoLib::GEOObjects and additionally emits signals when
 * data objects are modified. The GUI connects to these signals. Model instances
 * are created for every data object.
 */
class GEOModels : public QObject
{
	Q_OBJECT

public:
	GEOModels(QObject* parent = 0);
	~GEOModels();

	void setGEOObjects(GeoLib::GEOObjects* geo_objects);

	GeoTreeModel* getGeoModel() { return _geoModel; }
	StationTreeModel* getStationModel() { return _stationModel; }

public slots:
	/**
	 * Updates the tree model if the underlying data in GEOObjects has changed.
	 * Technically it removes the geometry from the tree model (but not from GeoObjects)
	 * and re-adds the (modified) geometry.
	 */
	void updateGeometry(const std::string &geo_name);

	/// Removes all parts (points, lines, surfaces) of the geometry with the given name.
	virtual void removeGeometry(std::string geo_name, GeoLib::GEOTYPE type);

	void addPointVec(std::string const& name);

	void removePointVec(std::string const& name);

	void addStationVec(std::string const& name);

	void removeStationVec(std::string const& name);

	void addPolylineVec(std::string const& name);

	void appendPolylineVec(std::string const& name);

	void removePolylineVec(std::string const& name);

	void addSurfaceVec(std::string const& name);

	/// @brief
	/// @param surfaces The surface vector.
	void appendSurfaceVec(std::string const& name);
	void removeSurfaceVec(std::string const& name);

	/// Adds the name 'new_name' for the geo-object specified by the parameters
	void addNameForElement(const std::string &geometry_name, const GeoLib::GEOTYPE object_type, std::size_t id, std::string new_name);

	/// Adds a generic name to all points that are part of the specified geo-object
	void addNameForObjectPoints(const std::string &geometry_name, const GeoLib::GEOTYPE object_type, const std::string &geo_object_name, const std::string &new_name);

	/// Calls all necessary functions to connect polyline-segments and update all views and windows.
	void connectPolylineSegments(const std::string &geoName,
	                             std::vector<std::size_t> indexlist,
	                             double proximity,
	                             std::string ply_name,
	                             bool closePly,
	                             bool triangulatePly);

protected:
	GeoTreeModel* _geoModel;
	StationTreeModel* _stationModel;

private:
	GeoLib::GEOObjects* _geo_objects;

signals:
	void geoDataAdded(GeoTreeModel*, std::string, GeoLib::GEOTYPE);
	void geoDataRemoved(GeoTreeModel*, std::string, GeoLib::GEOTYPE);

	void stationVectorAdded(StationTreeModel* model, std::string name);
	void stationVectorRemoved(StationTreeModel* model, std::string name);
};

class GEOModelsCallbacks final : public GeoLib::GEOObjects::Callbacks
{
public:
	explicit GEOModelsCallbacks(GEOModels& geo_models) : _geo_models(geo_models)
	{
	}

	void addPointVec(std::string const& name) override
	{
		_geo_models.addPointVec(name);
	}

	void removePointVec(std::string const& name) override
	{
		_geo_models.removePointVec(name);
	}

	void addStationVec(std::string const& name) override
	{
		_geo_models.addStationVec(name);
	};

	void removeStationVec(std::string const& name) override
	{
		_geo_models.removeStationVec(name);
	};

	void addPolylineVec(std::string const& name) override
	{
		_geo_models.addPolylineVec(name);
	};

	void appendPolylineVec(std::string const& name) override
	{
		_geo_models.appendPolylineVec(name);
	};

	void removePolylineVec(std::string const& name) override
	{
		_geo_models.removePolylineVec(name);
	};

	void addSurfaceVec(std::string const& name) override
	{
		_geo_models.addSurfaceVec(name);
	};

	void appendSurfaceVec(std::string const& name) override
	{
		_geo_models.appendSurfaceVec(name);
	};

	void removeSurfaceVec(std::string const& name) override
	{
		_geo_models.removeSurfaceVec(name);
	};

private:
	GEOModels& _geo_models;
};


#endif // GEOMODELS_H
