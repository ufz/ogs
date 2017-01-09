/**
 * \file
 * \author Thomas Fischer
 * \date   2012-02-27
 * \brief  Implementation of the OctTree class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstdint>

#include "MathLib/Point3d.h"
#include "MathLib/MathTools.h"

namespace GeoLib {

/// @tparam POINT point data type the OctTree will use
/// @tparam MAX_POINTS maximum number of pointers of POINT in a leaf
template <typename POINT, std::size_t MAX_POINTS>
class OctTree {
public:
    /// Create an OctTree object. The arguments ll and ur are used to compute a
    /// cube domain the OctTree will living in.
    /// @attention This cubic domain can not be resized during the life time of
    /// the OctTree.
    /// @param ll lower left front point, used for computation of cubic domain
    /// @param ur upper right back point, used for computation of cubic domain
    /// @param eps the euclidean distance as a threshold to make objects unique
    /// [default std::numeric_limits<double>::epsilon()]
    /// Adding a new item to an already "filled" OctTree node results in a
    /// split of the OctTree node. The smaller this number is the more leaves
    /// the OctTree will have, i.e. it needs more memory and more time to walk
    /// through the OctTree, but the search inside a leaf is fast. In
    /// contrast a big value results into a smaller number of OctTree leaves,
    /// the memory requirements for the OctTree may be lower but the search
    /// inside a OctTree leaf may be more expensive. The value should be
    /// choosen application dependend. [default 8]
    template <typename T>
    static OctTree<POINT, MAX_POINTS>* createOctTree(T ll, T ur,
        double eps = std::numeric_limits<double>::epsilon());

    /// Destroys the children of this node. @attention Does not destroy the
    /// pointers to the managed objects.
    virtual ~OctTree();

    /// This method adds the given point to the OctTree. If necessary,
    /// new OctTree nodes will be inserted deploying a split of the
    /// corresponding OctTree node.
    /// @param pnt the pointer to a point that should be inserted
    /// @param ret_pnt the pointer to a point in the OctTree. Three cases can
    /// occure:
    /// (1) ret_pnt is nullptr: the given point (pnt) is outside of the OctTree
    /// domain
    /// (2) ret_pnt is equal to pnt: the point is added to the OctTree
    /// (3) In case ret_pnt is neither equal to pnt nor equal to nullptr,
    /// another item within the eps distance is already in the OctTree and the
    /// pointer to this object is returned.
    /// @return If the point can be inserted the method returns true, else false.
    bool addPoint(POINT * pnt, POINT *& ret_pnt);

    /// range query - returns all points inside the range [min[0], max[0]) x
    /// [min[1], max[1]) x [min[2], max[2])
    template <typename T>
    void
    getPointsInRange(T const& min, T const& max, std::vector<POINT*> &pnts) const;

#ifndef NDEBUG
    MathLib::Point3d const& getLowerLeftCornerPoint() const { return _ll; }
    MathLib::Point3d const& getUpperRightCornerPoint() const { return _ur; }
    OctTree<POINT, MAX_POINTS> const* getChild(std::size_t i) const
    {
        return _children[i];
    }
    std::vector<POINT*> const& getPointVector() const { return _pnts; }
#endif

private:
    /// private constructor
    /// @param ll lower left point
    /// @param ur upper right point
    /// @param eps the euclidean distance as a threshold to make objects unique
    OctTree(MathLib::Point3d const& ll, MathLib::Point3d const& ur, double eps);

    enum class Quadrant : std::int8_t {
        NEL = 0, //!< north east lower
        NWL, //!< north west lower
        SWL, //!< south west lower
        SEL, //!< south east lower
        NEU, //!< south west upper
        NWU, //!< south west upper
        SWU, //!< south west upper
        SEU //!< south east upper
    };

    /// This method tries to add the given point to the OctTree. If necessary
    /// for adding the point, new nodes will be inserted into the OctTree.
    /// @param pnt, ret_pnt see documentation of addPoint()
    /// @return If the point can be inserted the method returns true, else false.
    bool addPoint_(POINT * pnt, POINT *& ret_pnt);

    /**
     * This method adds the given point to the OctTree. If necessary,
     * the OctTree will be extended.
     * @param pnt the point
     * @return If the point can be inserted the method returns true, else false.
     */
    bool addPointToChild(POINT* pnt);

    /// Creates the child nodes of this leaf and distribute the points stored
    /// in _pnts to the children.
    /// @param pnt the pointer to the points that is responsible for the split
    void splitNode(POINT * pnt);

    /// checks if the given point pnt is outside of the OctTree node.
    /// @param pnt the point that check is performed on
    /// @return true if the point is outside of the OctTree node.
    bool isOutside(POINT * pnt) const;

    /// children are sorted:
    ///   _children[0] is north east lower child
    ///   _children[1] is north west lower child
    ///   _children[2] is south west lower child
    ///   _children[3] is south east lower child
    ///   _children[4] is north east upper child
    ///   _children[5] is north west upper child
    ///   _children[6] is south west upper child
    ///   _children[7] is south east upper child
    std::array<OctTree<POINT, MAX_POINTS>*, 8> _children;
    /// lower left front face point of the cube
    MathLib::Point3d const _ll;
    /// upper right back face point of the cube
    MathLib::Point3d const _ur;
    /// vector of pointers to POINT objects
    std::vector<POINT*> _pnts;
    /// flag if this OctTree is a leaf
    bool _is_leaf;
    /// threshold for point uniqueness
    double const _eps;
};

} // end namespace GeoLib

#include "OctTree-impl.h"
