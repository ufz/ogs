/**
 * \file
 * \author Thomas Fischer
 * \date   2010-04-22
 * \brief  Definition of the AABB class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <bitset>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <tuple>
#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"
#include "MathLib/Point3d.h"
#include "MathLib/MathTools.h"

namespace GeoLib
{
/**
 * \brief Class AABB is an axis aligned bounding box around a given
 * set of geometric points of (template) type PNT_TYPE.
 *
 * Let \f$P = \{p_k \in \mathbb{R}^3, \ k=1, \dotsc, n\}\f$ a set of 3d points.
 * The bounding volume is described by its lower, left, front point \f$\ell\f$
 * and the upper, right, back point \f$u\f$, i.e. the coordinates of \f$\ell\f$
 * and \f$u\f$ are computed as follows \f$\ell_i = \min \limits_{p \in P}
 * p_i\f$, \f$u_i = \max \limits_{p \in P} p_i\f$, respectively. The bounding
 * box consists of half-open intervals \f$[\ell_x, u_x) \times [\ell_y, u_y)
 * \times [\ell_z, u_z).\f$ The bounding box is enlarged up to the next
 * available floating point number such that all input points are contained in
 * the bounding box.
 *
 */
class AABB
{
public:
    /**
     * construction of object, initialization the axis aligned bounding box
     * @tparam PNT_TYPE a point type supporting accessing the coordinates via
     * operator[]
     * */
    template <typename PNT_TYPE>
    AABB(std::vector<PNT_TYPE*> const& pnts, std::vector<std::size_t> const& ids)
    {
        assert(! ids.empty());
        init(pnts[ids[0]]);
        for (std::size_t i=1; i<ids.size(); ++i) {
            updateWithoutEnlarge(*(pnts[ids[i]]));
        }
        enlarge();
    }

    AABB(AABB const&) = default;

    /**
     * Construction of object using input iterators. In contrast to give a vector
     * this approach is more generic. You can use every (stl) container and
     * C arrays as input for constructing the object.
     * @attention{The constructor requires that std::distance(first, last) > 0.}
     * @param first the input iterator to the initial position in the sequence
     * @param last the input iterator to the final position in a container, i.e. [first, last).
     * @attention{The iterator last must be reachable from first.}
     */
    template <typename InputIterator>
    AABB(InputIterator first, InputIterator last)
    {
        if (std::distance(first,last) <= 0)
        {
            OGS_FATAL("AABB::AABB(InputIterator first, InputIterator last): first > last");
        }
        init(*first);
        InputIterator it(first);
        while (it != last) {
            updateWithoutEnlarge(*it);
            it++;
        }
        enlarge();
    }

    /// Checks if the bounding box has to be updated.
    /// @return true if AABB is updated.
    template <typename PNT_TYPE>
    bool update(PNT_TYPE const & p)
    {
        // First component of the pair signals if the minimum point is changed
        // Second component signals not only if the max point is changed.
        // Furthermore it is signaled what coordinate (0,1,2) is changed.
        std::pair<bool, std::bitset<3>> updated(false, 0);
        for (std::size_t k(0); k<3; k++) {
            // if the minimum point is updated pair.first==true
            if (p[k] < _min_pnt[k]) {
                _min_pnt[k] = p[k];
                updated.first = true;
            }
            // if the kth coordinate of the maximum point is updated
            // pair.second[k]==true
            if (p[k] >= _max_pnt[k]) {
                _max_pnt[k] = p[k];
                updated.second[k] = true;
            }
        }

        if (updated.second.any()) {
            enlarge(updated.second);
            return true;
        }
        return updated.first;
    }

    /**
     * check if point is in the axis aligned bounding box
     */
    template <typename T>
    bool containsPoint(T const & pnt, double eps) const
    {
        if (pnt[0] < _min_pnt[0] - eps || _max_pnt[0] + eps <= pnt[0])
            return false;
        if (pnt[1] < _min_pnt[1] - eps || _max_pnt[1] + eps <= pnt[1])
            return false;
        if (pnt[2] < _min_pnt[2] - eps || _max_pnt[2] + eps <= pnt[2])
            return false;
        return true;
    }

    template <typename T>
    bool containsPointXY(T const & pnt) const
    {
        if (pnt[0] < _min_pnt[0] || _max_pnt[0] <= pnt[0]) return false;
        if (pnt[1] < _min_pnt[1] || _max_pnt[1] <= pnt[1]) return false;
        return true;
    }

    /**
     * returns a point that coordinates are minimal for each dimension
     * for the given point set
     * @return a point
     */
    MathLib::Point3d const& getMinPoint() const { return _min_pnt; }

    /**
     * returns a point that coordinates are maximal for each dimension
     * within the given point set
     * @return a point
     */
    MathLib::Point3d const& getMaxPoint() const { return _max_pnt; }

    /**
     * Method checks if the given AABB object is contained within the
     * AABB represented by this object.
     * @param other_aabb the AABB to test with
     * @return true if the other AABB is contained in the AABB
     * represented by this object
     */
    bool containsAABB(AABB const& other_aabb) const
    {
        return containsPoint(other_aabb.getMinPoint(), 0) &&
               containsPoint(other_aabb.getMaxPoint(), 0);
    }

protected:
    MathLib::Point3d _min_pnt = MathLib::Point3d{std::array<double,3>{{
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()}}};
    MathLib::Point3d _max_pnt = MathLib::Point3d{std::array<double,3>{{
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest()}}};
private:
    /// Enlarge the bounding box the smallest possible amount (modifying the
    /// unit in the last place). Only the coordinates of the maximum point are
    /// changed such that the half-open property will be preserved.
    void enlarge(std::bitset<3> to_update = 7)
    {
        for (std::size_t k=0; k<3; ++k) {
            if (to_update[k]) {
                _max_pnt[k] = std::nextafter(_max_pnt[k],
                    std::numeric_limits<double>::max());
            }
        }
    }

    template <typename PNT_TYPE>
    void init(PNT_TYPE const & pnt)
    {
        _min_pnt[0] = _max_pnt[0] = pnt[0];
        _min_pnt[1] = _max_pnt[1] = pnt[1];
        _min_pnt[2] = _max_pnt[2] = pnt[2];
    }

    template <typename PNT_TYPE>
    void init(PNT_TYPE * const & pnt)
    {
        init(*pnt);
    }

    /// Private method that is used internally to update the min and max point
    /// of the bounding box using point \f$p\f$ without enlarging the bounding
    /// box. Using this method the bounding box of the initial point set is
    /// enlarged only once.
    /// @param p point that will possibly change the bounding box points
    template <typename PNT_TYPE>
    void  updateWithoutEnlarge(PNT_TYPE const & p)
    {
        for (std::size_t k(0); k<3; k++) {
            if (p[k] < _min_pnt[k]) {
                _min_pnt[k] = p[k];
            }
            if (p[k] >= _max_pnt[k]) {
                _max_pnt[k] = p[k];
            }
        }
    }

    template <typename PNT_TYPE>
    void updateWithoutEnlarge(PNT_TYPE * const & pnt)
    {
        updateWithoutEnlarge(*pnt);
    }

    template <typename PNT_TYPE>
    void update(PNT_TYPE const * pnt)
    {
        update(*pnt);
    }
};
} // end namespace
