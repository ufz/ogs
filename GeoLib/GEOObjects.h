/**
 * \file GEOObjects.h
 *
 *  Created on 2010-01-21 by Thomas Fischer / KR
 */

#ifndef GEOOBJECTS_H_
#define GEOOBJECTS_H_

#include <map>
#include <string>
#include <vector>

#include "Point.h"
#include "PointVec.h"
#include "Polyline.h"
#include "PolylineVec.h"
#include "Surface.h"
#include "SurfaceVec.h"

#include "Color.h"
#include "Station.h"

#include "GeoType.h"

namespace GeoLib
{
///
/**
 * \defgroup GeoLib This module consists of classes governing geometric objects
 * and related algorithms.
 */

/**
 * \brief Container class for geometric objects.
 *
 * This class contains all the methods necessary for the I/O of geometric objects.
 * Said objects are Points, polylines, Surfaces and Stations and they are stored in
 * vectors (arrays) which are identified by a unique name.
 * For a hierarchical definition, surfaces are bounded by polylines and polylines
 * are defined by points. Therefore, a vector of surfaces references a vector polylines
 * and a vector of polylines references a vector of points, respectively. For
 * identification purposes, all of these vectors have the same name, i.e. the polyline-
 * vector named "aaa" references a point vector "aaa". However, this name ("aaa") is
 * unique among all the vectors of the same class, i.e. there exists only one point-
 * vector with this name, etc.
 * Note: The fact that vectors are uniquely named and the same name is assigned to
 * related objects is automatically handled by this class.
 *
 * For each of these object-classes exists an "add", "remove" and "get"-method which
 * allows for loading/unloading as well as accessing the data, respectively.
 * E.g. for points these methods are "addPointVec(name)", "getPointVec(name)" and
 * "removePointVec(name)". For some objects, additional methods might exist if
 * necessary.
 */
class GEOObjects
{
public:
	/**
	 * Adds a vector of points with the given name to GEOObjects.
	 * @param points vector of pointers to points
	 * @param name the project name
	 * @param pnt_names vector of the names corresponding to the points
	 * @param eps relative tolerance value for testing of point uniqueness
	 */
	virtual void addPointVec(std::vector<Point*>* points,
	                         std::string &name,
	                         std::map<std::string, size_t>* pnt_names = NULL,
	                         double eps = sqrt(std::numeric_limits<double>::min()));

	/** copies the pointers to the points in the given vector to the PointVec of the provided name.
	 * The pointers are managed by the GEOObjects, i.e. GEOObjects will delete the Points at the
	 * end of its scope.
	 * \param points the vector with points
	 * \param name the name of the internal PointVec
	 * \param ids On return the vector holds the ids of the inserted points within the internal vector.
	 * In case the ids are not needed it is possible to give a NULL pointer (default value).
	 * \return true if the points are appended, false if the a PointVec with the
	 * corresponding name does not exist
	 * */
	virtual bool appendPointVec(const std::vector<Point*> &points,
	                            std::string const &name, std::vector<size_t>* ids = NULL);

	/**
	 * Method appends the point the the PointVec object with the name name. The PointVec
	 * object takes care about deleting the point. If the point already exists within the
	 * PointVec object this method will delete it.
	 * @param point (input) the point (exact the pointer to the point) that should be added
	 * @param name (input) the name of the geometry the point should be added
	 * @param id (output) the id of the point within the PointVec object will be set
	 * @return true, if the point could be inserted (i.e. was not already contained in the vector)
	 * else false (the user have to delete the point itself)
	 */
	bool appendPoint(Point* point, std::string const &name, size_t& id);

	/**
	 * Returns the point vector with the given name.
	 */
	const std::vector<Point*>* getPointVec(const std::string &name) const;

	/**
	 * search and returns the PointVec object with the given name.
	 * @param name the name of the PointVec object
	 * @return the PointVec object stored in GEOObjects
	 */
	const PointVec* getPointVecObj(const std::string &name) const;

	/// Returns a pointer to a PointVec object for the given name.
	PointVec* getPointVecObj(const std::string &name)
	{
		return const_cast<PointVec*>(static_cast<const GEOObjects&>(*this).
		                                getPointVecObj(name));
	}

	/** If there exists no dependencies the point vector with the given
	 * name from GEOObjects will be removed and the method returns true,
	 * else the return value is false.
	 */
	virtual bool removePointVec(const std::string &name);

	/// Adds a vector of stations with the given name and colour to GEOObjects.
	virtual void addStationVec(std::vector<Point*>* stations, std::string &name);

	/// Filters a list of stations with the given name based on the criteria in PropertyBounds.
	/// (See property system in Station class for more information.)
	std::vector<Point*>* filterStationVec(const std::string &name,
	                                      const std::vector<PropertyBounds> &bounds);
	/// Returns the station vector with the given name.
	const std::vector<Point*>* getStationVec(const std::string &name) const;

	/// Removes the station vector with the given name from GEOObjects
	virtual bool removeStationVec(const std::string &name)
	{
		return removePointVec(name);
	}

	/**
	 * Adds a vector of polylines with the given name to GEOObjects.
	 * @param lines The lines vector.
	 * @param name The given name.
	 * @param ply_names vector of the names corresponding to the polylines
	 */
	virtual void addPolylineVec(std::vector<Polyline*>* lines,
	                            const std::string &name,
	                            std::map<std::string,size_t>* ply_names = NULL);

	/** copies the pointers to the polylines in the vector to the PolylineVec with provided name.
	 * the pointers are managed by the GEOObjects, i.e. GEOObjects will delete the Polylines at the
	 * end of its scope
	 * \param polylines the vector with polylines
	 * \param name the name of the internal PolylineVec
	 * \return true if the polylines are appended, false if the PolylineVec with the
	 * corresponding name does not exist
	 * */
	virtual bool appendPolylineVec(const std::vector<Polyline*> &polylines,
	                               const std::string &name);

	/**
	 * Returns the polyline vector with the given name.
	 * */
	const std::vector<Polyline*>* getPolylineVec(const std::string &name) const;

	/**
	 * Returns a pointer to a PolylineVec object for the given name as a const.
	 * @param name the name of the vector of polylines
	 * @return PolylineVec object
	 */
	const PolylineVec* getPolylineVecObj(const std::string &name) const;

	/// Returns a pointer to a PolylineVec object for the given name.
	PolylineVec* getPolylineVecObj(const std::string &name)
	{
		return const_cast<PolylineVec*>(static_cast<const GEOObjects&>(*this).
		                                getPolylineVecObj(name));
	}

	/**
	 * If no Surfaces depends on the vector of Polylines with the given
	 * name it will be removed and the method returns true,
	 * else the return value is false.
	 */
	virtual bool removePolylineVec(const std::string &name);

	/** Adds a vector of surfaces with the given name to GEOObjects. */
	virtual void addSurfaceVec(std::vector<Surface*>* surfaces,
	                           const std::string &name,
	                           std::map<std::string, size_t>* sfc_names = NULL);

	/**
	 * Copies the surfaces in the vector to the SurfaceVec with the given name.
	 * \param surfaces the vector with surfaces
	 * \param name the name of the internal PolylineVec
	 * \return true if the surfaces are appended, false if the SurfaceVec with the
	 * corresponding name does not exist
	 * */
	virtual bool appendSurfaceVec(const std::vector<Surface*> &surfaces,
	                              const std::string &name);

	/// Returns the surface vector with the given name as a const.
	const std::vector<Surface*>* getSurfaceVec(const std::string &name) const;

	/// Returns the surface vector with the given name.
	SurfaceVec* getSurfaceVecObj(const std::string &name)
	{
		return const_cast<SurfaceVec*>(static_cast<const GEOObjects&>(*this).
		                               getSurfaceVecObj(name));
	}

	/** removes the vector of Surfaces with the given name */
	virtual bool removeSurfaceVec(const std::string &name);
	/**
	 * Returns a pointer to a SurfaceVec object for the given name. The class
	 * SurfaceVec stores the relation between surfaces and the names of the surfaces.
	 * @param name the name of the vector of surfaces (the project name)
	 * @return SurfaceVec object
	 */
	const SurfaceVec* getSurfaceVecObj(const std::string &name) const;

	/// Returns -1 if no geometry of the given name exists or die index of the geometry in _pnt_vecs otherwise
	int exists(const std::string &geometry_name) const;

	/// Returns the names of all geometry vectors.
	void getGeometryNames (std::vector<std::string>& names) const;

	const std::string getElementNameByID(const std::string &geometry_name, GeoLib::GEOTYPE type, size_t id) const;

	/// Returns the names of all station vectors.
	void getStationVectorNames(std::vector<std::string>& names) const;

	/**
	 * Method mergeGeometries merges the geometries that are given by the names in the vector.
	 * Stations points are not included in the resulting merged geometry.
	 * @param names the names of the geometries that are to be merged
	 * @param merged_geo_name the name of the resulting geometry
	 */
	void mergeGeometries(std::vector<std::string> const & names, std::string &merged_geo_name);

	/// Returns the geo object for a geometric item of the given name and type for the associated geometry.
	const GeoLib::GeoObject* getGEOObject(const std::string &geo_name,
	                                      GeoLib::GEOTYPE type,
	                                      const std::string &obj_name) const;

	/** constructor */
	GEOObjects();
	/** destructor */
	virtual ~GEOObjects();

protected:
	/**
	 * Determines if the given name is unique among all the names in point vectors and creates a
	 * new name if this is not the case. The new name is then simply "name + x", where x>1 is
	 * the smallest number that creates a unique name (i.e. "name-2", "name-3", etc.)
	 * \param name Original name of the list, this name might be changed within this method if necessary.
	 * \return true if the name was unique, false if a new name has been generated
	 */
	bool isUniquePointVecName(std::string &name);

	/// Checks if the point vector with the given name is referenced in a polyline- or surface vector.
	bool isPntVecUsed (const std::string &name) const;

	/**
	 * vector manages pointers to PointVec objects
	 */
	std::vector<PointVec*> _pnt_vecs;

	/** vector manages pointers to PolylineVec objects */
	std::vector<PolylineVec*> _ply_vecs;
	/** vector manages pointers to SurfaceVec objects */
	std::vector<SurfaceVec*> _sfc_vecs;
};
} // end namespace

#endif /* GEOOBJECTS_H_ */
