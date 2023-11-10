/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-14
 * \brief  Definition of the PolyLine class.
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cmath>
#include <vector>

// GeoLib
#include "GeoObject.h"
#include "LineSegment.h"
#include "Point.h"

// MathLib
#include "MathLib/Point3d.h"

namespace GeoLib
{
class PointVec;

/**
 * \ingroup GeoLib
 *
 * \brief Class Polyline consists mainly of a reference to a point vector and
 * a vector that stores the indices in the point vector.
 * A polyline consists of at least one line segment. The polyline is specified
 * by the points of the line segments. The class Polyline stores ids of pointers
 * to the points in the _ply_pnt_ids vector.
 * */
class Polyline : public GeoObject
{
public:
    class SegmentIterator final
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = LineSegment;
        using difference_type = std::ptrdiff_t;
        using pointer = LineSegment*;
        using reference = LineSegment&;

        explicit SegmentIterator(Polyline const& polyline,
                                 std::size_t segment_number);

        SegmentIterator(SegmentIterator const& src);

        SegmentIterator() = delete;
        ~SegmentIterator() = default;

        SegmentIterator& operator=(SegmentIterator const& rhs);

        std::size_t getSegmentNumber() const;

        SegmentIterator& operator++();

        LineSegment operator*() const;

        LineSegment operator*();

        bool operator==(SegmentIterator const& other) const;

        bool operator!=(SegmentIterator const& other) const;

        SegmentIterator& operator+=(
            std::vector<GeoLib::Point>::difference_type n);

        SegmentIterator operator+(
            std::vector<GeoLib::Point>::difference_type n);

        SegmentIterator& operator-=(
            std::vector<GeoLib::Point>::difference_type n);

        SegmentIterator operator-(
            std::vector<GeoLib::Point>::difference_type n);

    private:
        GeoLib::Polyline const* _polyline;
        std::vector<GeoLib::Point*>::size_type _segment_number;
    };

    friend class Polygon;
    /** constructor
     * \param pnt_vec a reference to the point vector
     */
    explicit Polyline(const std::vector<Point*>& pnt_vec);
    /**
     * Copy constructor
     * \param ply Polyline
     */
    Polyline(const Polyline& ply);
    Polyline& operator=(Polyline const& other) = delete;

    ~Polyline() override = default;

    /// return a geometry type
    GEOTYPE getGeoType() const override { return GEOTYPE::POLYLINE; }

    /**
     * Adds an id of a point at the end of the polyline if and only if the
     * resulting segment won't be empty. The id have to be inside the (internal)
     * \c _ply_pnts vector the polyline is based on.
     * \return If the point could be added the return value is \c true. If the
     * addition of the point would result in empty line segment \c false is
     * returned.
     */
    virtual bool addPoint(std::size_t pnt_id);

    /**
     * Method inserts a new point (that have to be inside the _ply_pnts vector)
     * at the given position in the polyline if and only if the resulting
     * segments won't be empty.
     * \param pos the position in the polyline, pos have to be a value into the
     * interval [0, number of points)
     * \param pnt_id the id of the new point in the vector of points the
     * polyline is based on
     * \return true if the point could be inserted, else false (if empty line
     * segments would be created).
     */
    virtual bool insertPoint(std::size_t pos, std::size_t pnt_id);

    /**
     * Method removes a point from the polyline. The connecting line segments
     * will be removed and the length of the polyline will be changed.
     * \param pos a valid position within the polyline
     */
    void removePoint(std::size_t pos);

    /**
     * Closes a polyline by adding a line segment that connects start- and
     * end-point.
     */
    void closePolyline();

    /// Constructs one polyline from a vector of connected polylines.
    /// All polylines in this vector need to reference the same point vector.
    static Polyline* constructPolylineFromSegments(
        const std::vector<Polyline*>& ply_vec, double prox = 0.0);

    /**
     * returns the number of points,
     * the number of segments is about one smaller
     * */
    std::size_t getNumberOfPoints() const;

    std::size_t getNumberOfSegments() const;

    /** returns true if the polyline is closed */
    bool isClosed() const;

    /** returns true if the polyline is coplanar */
    bool isCoplanar() const;

    /**
     * Method tests if the given id of a point (within the vector of points the
     * polyline is based on) is inside the polyline.
     * \param pnt_id the id of the point
     * \return true if the point is part of the polyline, else false
     */
    bool isPointIDInPolyline(std::size_t pnt_id) const;

    /**
     * returns the index of the i-th polyline point
     * in the point vector
     */
    std::size_t getPointID(std::size_t const i) const;

    /**
     * Changes a point index for one point in a line
     * \param idx Index of point in line
     * \param id ID of point in PointVec object
     */
    void setPointID(std::size_t idx, std::size_t id);

    /**
     * \brief returns the i-th point contained in the polyline
     * */
    const Point* getPoint(std::size_t i) const;

    SegmentIterator begin() const { return SegmentIterator(*this, 0); }

    SegmentIterator end() const
    {
        return SegmentIterator(*this, getNumberOfSegments());
    }

    std::vector<Point*> const& getPointsVec() const;

    /**
     * returns the distance along the polyline from the beginning of the
     * polyline
     * \param pnt the point on the polyline
     * \param epsilon_radius the epsilon
     * \return the distance along the polyline between the given point and the
     * beginning of the polyline. If the given point is not on the
     * polyine, negative value is returned.
     */
    double getDistanceAlongPolyline(const MathLib::Point3d& pnt,
                                    const double epsilon_radius) const;

protected:
    /** a reference to the vector of pointers to the geometric points */
    const std::vector<Point*>& _ply_pnts;

    void reverseOrientation();

    std::vector<std::size_t> const& getPolylinePointIDs() const
    {
        return _ply_pnt_ids;
    }

private:
    /** position of pointers to the geometric points */
    std::vector<std::size_t> _ply_pnt_ids;

    LineSegment getSegment(std::size_t i) const;
};

bool containsEdge(const Polyline& ply, std::size_t id0, std::size_t id1);

/// Resets the point IDs of the polyline corresponding to the mapping.
void resetPointIDs(Polyline& polyline, std::vector<std::size_t> const& mapping);

/// Resets the point IDs of the polyline corresponding to the mapping.
void markUsedPoints(Polyline const& polyline, std::vector<bool>& used_points);

/**
 * comparison operator
 * \param lhs first polyline
 * \param rhs second polyline
 * \return true, if the polylines consists of the same sequence of line segments
 */
bool operator==(Polyline const& lhs, Polyline const& rhs);

bool pointsAreIdentical(const std::vector<Point*>& pnt_vec, std::size_t i,
                        std::size_t j, double prox);

/// Create a polyline from given point ids.
std::unique_ptr<Polyline> createPolyline(GeoLib::PointVec const& points_vec,
                                         std::vector<std::size_t>&& point_ids);

}  // namespace GeoLib
