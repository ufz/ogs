/**
 * \file
 * \author Thomas Fischer / Karsten Rink
 * \date   2010-02-02
 * \brief  Definition of the PointVec class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "AABB.h"
#include "OctTree.h"
#include "Point.h"
#include "TemplateVec.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 *
 * \brief This class manages pointers to Points in a std::vector along with a
 * name. It also handles the deletion of points. Furthermore, typically,
 * PointVec objects are managed by GEOObjects using the instance name for
 * identification. For this reason PointVec must have a unique name.
 * */
class PointVec final : public TemplateVec<Point>
{
public:
    /// Signals if the vector contains object of type Point or Station
    enum class PointType
    {
        POINT = 0,
        STATION = 1
    };

    /**
     * Constructor initializes the name of the PointVec object,
     * and moves the points vector into the internal data structure.
     * @param name the name of the point set
     * @param points R-value reference to a vector of pointers to
     * GeoLib::Points.
     * @attention{The PointVec object takes the ownership of the GeoLib::Points
     * in the points vector, i.e. it delete the points during destruction!}
     * @param name_id_map A std::map that stores the relation name to point.
     * @param type the type of the point, \sa enum PointType
     * @param rel_eps This is a relative error tolerance value for the test of
     * identical points. The size of the axis aligned bounding box multiplied
     * with the value of rel_eps gives the real tolerance \f$tol\f$. Two points
     * \f$p_0, p_1 \f$ are identical iff \f$|p_1 - p_0| \le tol.\f$
     */
    PointVec(std::string const& name, std::vector<Point*>&& points,
             NameIdMap&& name_id_map,
             PointType const type = PointVec::PointType::POINT,
             double const rel_eps = std::numeric_limits<double>::epsilon());
    /**
     * This constructor is similar to the above constructor with the
     * exception that the NameIdMap must not be passed.
     */
    PointVec(std::string const& name, std::vector<Point*>&& points,
             PointType const type = PointVec::PointType::POINT,
             double const rel_eps = std::numeric_limits<double>::epsilon());

    /**
     * Method adds a Point to the (internal) standard vector and takes the
     * ownership. If the given point is already included in the vector, the
     * point will be destroyed and the id of the existing point will be
     * returned.
     * @param pnt the pointer to the Point
     * @return the id of the point within the internal vector
     */
    std::size_t push_back(Point* pnt);

    /**
     * push_back adds new elements at the end of the vector _data_vec.
     * @param pnt a pointer to the point, PointVec takes ownership of the point
     * @param name the name of the point
     */
    void push_back(Point* pnt, std::string const* const name) override;

    /**
     * get the type of Point, this can be either POINT or STATION
     */
    PointType getType() const { return _type; }

    const std::vector<std::size_t>& getIDMap() const { return _pnt_id_map; }

    const GeoLib::AABB& getAABB() const;

    std::string const& getItemNameByID(std::size_t id) const;

    void setNameForElement(std::size_t id, std::string const& name) override;

    /// Resets the internal data structures, i.e., the axis aligned bounding
    /// box, the relative epsilon for the equality tests and the oct tree data
    /// structure.
    /// \note This method have to be called if the coordinates of a point stored
    /// by the PointVec is modified from outside.
    void resetInternalDataStructures();

private:
    /**
     * After the point set is modified (for example by makePntsUnique()) the
     * mapping has to be corrected.
     */
    void correctNameIDMapping();

    PointVec(PointVec&&) = delete;
    PointVec(const PointVec&) = delete;
    PointVec() = delete;
    PointVec& operator=(const PointVec& rhs) = delete;
    PointVec& operator=(PointVec&& rhs) = delete;

    /**
     * Inserts the instance of the Point into internal data structures
     * (@see TemplateVec::_data_vec, _pnt_id_map) if and only if there
     * does not exist a point with the same coordinates (up to
     * std::numeric_limits<double>::epsilon()). In case there exists
     * already a point with the same coordinates the given pnt-object
     * will be deleted!
     * @param pnt  Pointer to GeooLib::Point instance
     * @return either the new id or the id of the existing point with
     * the same coordinates
     */
    std::size_t uniqueInsert(Point* pnt);

    /** the type of the point (\sa enum PointType) */
    PointType _type;

    /**
     * permutation of the geometric elements according
     * to their lexicographical order
     */
    std::vector<std::size_t> _pnt_id_map;

    /// The reverse map to the name to id map, for fast lookup of the name to a
    /// given point id.
    std::vector<std::string> _id_to_name_map;

    AABB _aabb;
    double _rel_eps;
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 16>> _oct_tree;
};
}  // namespace GeoLib
