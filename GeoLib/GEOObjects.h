/**
 * \author Thomas Fischer / Karsten Rink
 * \date   2010-01-21
 * \brief  Definition of the GEOObjects class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <functional>
#include <string>
#include <vector>

#include "Point.h"
#include "PointVec.h"
#include "Polyline.h"
#include "PolylineVec.h"
#include "Surface.h"
#include "SurfaceVec.h"

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
class GEOObjects final
{
public:
    struct Callbacks
    {
        virtual void addPointVec(std::string const&){}
        virtual void removePointVec(std::string const&){}
        virtual void addStationVec(std::string const&){}
        virtual void removeStationVec(std::string const&){}
        virtual void addPolylineVec(std::string const&){}
        virtual void appendPolylineVec(std::string const&){}
        virtual void removePolylineVec(std::string const&){}
        virtual void addSurfaceVec(std::string const&){}
        virtual void appendSurfaceVec(std::string const&){}
        virtual void removeSurfaceVec(std::string const&){}
        virtual void renameGeometry(std::string const&, std::string const&) {}
        virtual ~Callbacks() = default;
    };

public:
    /**
     * Adds a vector of points with the given name to GEOObjects.
     * @param points vector of pointers to points
     * @param name the project name
     * @param pnt_names vector of the names corresponding to the points
     * @param eps relative tolerance value for testing of point uniqueness
     */
    void addPointVec(
        std::unique_ptr<std::vector<Point*>> points,
        std::string& name,
        std::unique_ptr<std::map<std::string, std::size_t>> pnt_names = nullptr,
        double eps = sqrt(std::numeric_limits<double>::epsilon()));

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
    bool removePointVec(const std::string &name);

    /// Adds a vector of stations with the given name and colour to GEOObjects.
    void addStationVec(
        std::unique_ptr<std::vector<Point *>> stations, std::string &name);

    /// Returns the station vector with the given name.
    const std::vector<GeoLib::Point*>* getStationVec(
        const std::string& name) const;

    /// Removes the station vector with the given name from GEOObjects
    bool removeStationVec(const std::string &name)
    {
        _callbacks->removeStationVec(name);
        return removePointVec(name);
    }

    /**
     * Adds a vector of polylines with the given name to GEOObjects.
     * @param lines The lines vector.
     * @param name The geometry to which the given Polyline objects should be added.
     * @param ply_names map of names and ids that are corresponding to the polylines
     */
    void addPolylineVec(std::unique_ptr<std::vector<Polyline*>> lines,
                        const std::string& name,
                        std::unique_ptr<std::map<std::string, std::size_t>>
                            ply_names = nullptr);

    /** copies the pointers to the polylines in the vector to the PolylineVec with provided name.
     * the pointers are managed by the GEOObjects, i.e. GEOObjects will delete the Polylines at the
     * end of its scope
     * \param polylines the vector with polylines
     * \param name the name of the internal PolylineVec
     * \return true if the polylines are appended, false if the PolylineVec with the
     * corresponding name does not exist
     * */
    bool appendPolylineVec(const std::vector<Polyline*> &polylines,
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
    bool removePolylineVec(const std::string &name);

    /** Adds a vector of surfaces with the given name to GEOObjects. */
    void addSurfaceVec(std::unique_ptr<std::vector<Surface*>> surfaces,
                       const std::string& name,
                       std::unique_ptr<std::map<std::string, std::size_t>>
                           sfc_names = nullptr);

    /**
     * Copies the surfaces in the vector to the SurfaceVec with the given name.
     * \param surfaces the vector with surfaces
     * \param name the name of the internal PolylineVec
     * \return true if the surfaces are appended, false if the SurfaceVec with the
     * corresponding name does not exist and the surfaces are added.
     * */
    bool appendSurfaceVec(const std::vector<Surface*>& surfaces,
                          const std::string& name);

    /// Returns the surface vector with the given name as a const.
    const std::vector<Surface*>* getSurfaceVec(const std::string &name) const;

    /// Returns the surface vector with the given name.
    SurfaceVec* getSurfaceVecObj(const std::string &name)
    {
        return const_cast<SurfaceVec*>(static_cast<const GEOObjects&>(*this).
                                       getSurfaceVecObj(name));
    }

    /** removes the vector of Surfaces with the given name */
    bool removeSurfaceVec(const std::string &name);
    /**
     * Returns a pointer to a SurfaceVec object for the given name. The class
     * SurfaceVec stores the relation between surfaces and the names of the surfaces.
     * @param name the name of the vector of surfaces (the project name)
     * @return SurfaceVec object
     */
    const SurfaceVec* getSurfaceVecObj(const std::string &name) const;

    /// Returns the names of all geometry vectors.
    void getGeometryNames (std::vector<std::string>& names) const;

    const std::string getElementNameByID(const std::string &geometry_name, GeoLib::GEOTYPE type, std::size_t id) const;

    /// Returns the names of all station vectors.
    void getStationVectorNames(std::vector<std::string>& names) const;

    /**
     * Determines if the given name is unique among all the names in point vectors and creates a
     * new name if this is not the case. The new name is then simply "name + x", where x>1 is
     * the smallest number that creates a unique name (i.e. "name-2", "name-3", etc.)
     * \param name Original name of the list, this name might be changed within this method if necessary.
     * \return true if the name was unique, false if a new name has been generated
     */
    bool isUniquePointVecName(std::string &name);

    /**
     * Method mergeGeometries merges the geometries that are given by the names in the vector.
     * Stations points are not included in the resulting merged geometry.
     * @param names the names of the geometries that are to be merged
     * @param merged_geo_name the name of the resulting geometry
     * @return 0 if success, 1 if no point-list is found for at least one of the geometries and 2 if the mergelist only contains less than two geometry
     */
    int mergeGeometries(std::vector<std::string> const & names, std::string &merged_geo_name);

    /// Renames an existing geometry, i.e. renames the internal PointVec,
    /// PolylineVec and the SurfaceVec objects from \c old_name to \c new_name.
    /// If no such PointVec, PolylineVec and SurfaceVec objects exist nothing
    /// will happen.
    void renameGeometry(std::string const& old_name,
                        std::string const& new_name);

    /// Returns the geo object for a geometric item of the given name and type for the associated geometry.
    const GeoLib::GeoObject* getGeoObject(const std::string &geo_name,
                                          GeoLib::GEOTYPE type,
                                          const std::string &geo_obj_name) const;

    /// Return named (by the tuple geo_name and geo_obj_name) geo object.
    // It is required that a tuple is a unique key for a geometric object!
    // If there is another geo object with same name one of them is returned.
    // In theory different types of geometric objects can have the same name.
    // For instance it is possible that a point object and a polyline object
    // share the same name. If there exists several objects sharing the same
    // name the first object found will be returned.
    // @param geo_name name of geometry
    // @param geo_obj_name name of the geo object
    GeoLib::GeoObject const* getGeoObject(const std::string &geo_name,
                                          const std::string &geo_obj_name) const;

    /** constructor */
    GEOObjects();
    /** destructor */
    ~GEOObjects();

    /// Returns std::numeric_limits<std::size_t>::max() if no geometry of the
    /// given name exists or the index of the geometry in _pnt_vecs otherwise
    std::size_t exists(const std::string &geometry_name) const;

    /// Checks if the point vector with the given name is referenced in a polyline- or surface vector.
    bool isPntVecUsed (const std::string &name) const;

    /// Read access to points w/o using a name.
    std::vector<PointVec*> const& getPoints() const { return _pnt_vecs; }
    /// Read access to polylines w/o using a name.
    std::vector<PolylineVec*> const& getPolylines() const { return _ply_vecs; }
    /// Read access to surfaces w/o using a name.
    std::vector<SurfaceVec*> const& getSurfaces() const { return _sfc_vecs; }

    /**
     * vector manages pointers to PointVec objects
     */
    std::vector<PointVec*> _pnt_vecs;

    /** vector manages pointers to PolylineVec objects */
    std::vector<PolylineVec*> _ply_vecs;
    /** vector manages pointers to SurfaceVec objects */
    std::vector<SurfaceVec*> _sfc_vecs;

    std::unique_ptr<Callbacks> _callbacks{new Callbacks};

    std::function<void(std::string const&)> addPolylineVecCallback =
        [](std::string const&)
    {
    };

    std::function<void(std::string const&)> appendPolylineVecCallback =
        [](std::string const&)
    {
    };

    std::function<void(std::string const&)> removePolylineVecCallback =
        [](std::string const&)
    {
    };

    std::function<void(std::string const&)> addSurfaceVecCallback =
        [](std::string const&)
    {
    };

    std::function<void(std::string const&)> appendSurfaceVecCallback =
        [](std::string const&)
    {
    };

    std::function<void(std::string const&)> removeSurfaceVecCallback =
        [](std::string const&)
    {
    };

private:
    /**
     * Method merges points from different geometries into one geometry. This
     * is a helper method for GEOObjects::mergeGeometries().
     * @param geo_names The vector of names of the geometries to merge.
     * @param merged_geo_name The (new) name of the geometry resulting from
     * merging.
     * @param pnt_offsets offsets in the merged vector storing the points
     * @return true, if merging the points succeeded, else false
     */
    bool mergePoints(std::vector<std::string> const & geo_names, std::string & merged_geo_name,
            std::vector<std::size_t> &pnt_offsets);

    /**
     * Method merges GeoLib::Polylines from different geometries into one
     * geometry. There isn't a check if the polyline is unique within the
     * given data sets. If there are two polylines with the same name, the
     * second occurrence will lost their name. This is a helper for
     * GEOObjects::mergeGeometries() and should be used only after
     * GEOObjects::mergePoints() is executed.
     * @param geo_names The vector of names of the geometries to merge.
     * @param merged_geo_name The (new) name of the geometry resulting from
     * merging.
     * @param pnt_offsets offsets in the merged vector storing the points.
     */
    void mergePolylines(std::vector<std::string> const & geo_names, std::string & merged_geo_name,
            std::vector<std::size_t> const& pnt_offsets);

    /**
     * Method merges GeoLib::Surfaces from different geometries into one
     * geometry. There isn't a check if the GeoLib::Surface is unique within
     * the given data sets. This is a helper for GEOObjects::mergeGeometries()
     * and should be used only after GEOObjects::mergePoints() is executed.
     * @param geo_names The vector of names of the geometries to merge.
     * @param merged_geo_name The (new) name of the geometry resulting from
     * merging.
     * @param pnt_offsets offsets in the merged vector storing the points.
     */
    void mergeSurfaces(std::vector<std::string> const & geo_names,
            std::string & merged_geo_name, std::vector<std::size_t> const& pnt_offsets);
};
} // end namespace
