/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <memory>

#include "GeoObject.h"
#include "Point.h"
#include "AABB.h"

namespace GeoLib {

class Polyline;

class Triangle;
class SurfaceGrid;

/**
 * \ingroup GeoLib
 *
 * \brief A Surface is represented by Triangles. It consists of a reference
 * to a vector of (pointers to) points (_sfc_pnts) and a vector that stores
 * the Triangles consisting of points from _sfc_pnts.
 * */
class Surface final : public GeoObject
{
public:
    explicit Surface(const std::vector<Point*>& pnt_vec);
    Surface(Surface const& src);
    ~Surface() override;

    Surface(Surface && src) = delete;
    Surface& operator=(Surface const& src) = delete;
    Surface& operator=(Surface && src) = delete;

    /// return a geometry type
    GEOTYPE getGeoType() const override { return GEOTYPE::SURFACE; }
    /**
     * adds three indices describing a triangle and updates the bounding box
     * */
    void addTriangle(std::size_t pnt_a, std::size_t pnt_b, std::size_t pnt_c);

    /**
     * returns the number of triangles describing the Surface
     * */
    std::size_t getNumberOfTriangles() const;

    /** \brief const access operator for the access to the i-th Triangle of the
     * surface.
    */
    const Triangle* operator[](std::size_t i) const;

    /**
     * is the given point in the bounding volume of the surface
     */
    bool isPntInBoundingVolume(MathLib::Point3d const& pnt, double eps) const;

    /**
     * is the given point pnt located in the surface
     * @param pnt the point
     * @param eps geometric tolerance for the search
     * @return true if the point is contained in the surface
     */
    bool isPntInSfc(MathLib::Point3d const& pnt, double eps) const;

    /**
     * find a triangle in which the given point is located
     * @param pnt the point
     * @return a pointer to a triangle. nullptr is returned if the point is not
     * contained in the surface
     */
    const Triangle* findTriangle(MathLib::Point3d const& pnt) const;

    const std::vector<Point*>* getPointVec() const { return &_sfc_pnts; }
    /**
     * method allows access to the internal axis aligned bounding box
     * @return axis aligned bounding box
     */
    AABB const& getAABB() const { return *_bounding_volume; }
protected:
    /** a vector of pointers to Points */
    const std::vector<Point*>& _sfc_pnts;
    /** position of pointers to the geometric points */
    std::vector<Triangle*> _sfc_triangles;
    /** bounding volume is an axis aligned bounding box */
    std::unique_ptr<AABB> _bounding_volume;
    /// The surface grid is a helper data structure to accelerate the point
    /// search. The method addTriangle() invalidates/resets the surface grid.
    /// A valid surface grid is created in case the const method isPntInSfc() is
    /// called and a valid surface grid is not existing.
    mutable std::unique_ptr<SurfaceGrid> _surface_grid;
};
}
