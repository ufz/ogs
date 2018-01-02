/**
 * \file
 * \author Thomas Fischer
 * \date   2010-11-09
 * \brief  Definition of the QuadTree class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>

#include <logog/include/logog.hpp>
#include <utility>

namespace GeoLib
{
/**
 * A quadtree is a rooted tree in which every internal
 * node has four children. Every node corresponds to a square.
 * (see Computational Geometry - Algorithms and Applications [Mark de Berg,
 * Otfried Cheong, Marc van Kreveld, Mark Overmars] - Chapter 14)
 *
 * One can instantiate the class template with a point type and a value for
 * the maximal number of points per node. The point-type have to provide
 * the access to its coordinates via operator[] and for debugging
 * purposes operator<<)
 */
template <typename POINT> class QuadTree
{
public:
    enum class Quadrant {
        NE = 0, //!< north east
        NW, //!< north west
        SW, //!< south west
        SE //!< south east
    };
    /**
     * This is the constructor for class QuadTree. It takes two points
     * (lower left and the upper right points).
     * @param ll lower left point of the square
     * @param ur upper right point of the square
     * @param max_points_per_leaf maximum number of points per leaf
     */
    QuadTree(POINT ll, POINT ur, std::size_t max_points_per_leaf)
        : _father(nullptr),
          _ll(std::move(ll)),
          _ur(std::move(ur)),
          _depth(0),
          _is_leaf(true),
          _max_points_per_leaf(max_points_per_leaf)
    {
        assert (_max_points_per_leaf > 0);

        // init children
        for (auto& child : _children)
            child = nullptr;

        if ((_ur[0] - _ll[0]) > (_ur[1] - _ll[1]))
            _ur[1] = _ll[1] + _ur[0] - _ll[0];
        else
            _ur[0] = _ll[0] + _ur[1] - _ll[1];

        DBUG("QuadTree(): lower left: (%f,%f,%f), upper right: (%f,%f,%f), depth: %d", _ll[0], _ll[1], _ll[2], _ur[0], _ur[1], _ur[2], _depth);
    }

    /**
     * destructor
     */
    ~QuadTree()
    {
        for (auto& child : _children)
        {
            delete child;
        }
    }

    /**
     * This method adds the given point to the quadtree. If necessary,
     * the quadtree will be extended.
     * @param pnt the point
     * @return If the point can be inserted the method returns true, else false.
     */
    bool addPoint (POINT const* pnt)
    {
        if ((*pnt)[0] < _ll[0]) return false;
        if ((*pnt)[0] >= _ur[0]) return false;
        if ((*pnt)[1] < _ll[1]) return false;
        if ((*pnt)[1] >= _ur[1]) return false;

        if (!_is_leaf) {
            for (auto& child : _children)
            {
                if (child->addPoint(pnt))
                    return true;
            }
            return false;
        }

        // check if point is already in quadtree
        bool pnt_in_quadtree (false);
        for (std::size_t k(0); k < _pnts.size() && !pnt_in_quadtree; k++) {
            const double v0((*(_pnts[k]))[0] - (*pnt)[0]);
            const double v1((*(_pnts[k]))[1] - (*pnt)[1]);
            const double sqr_dist (v0*v0 + v1*v1);
            if (sqr_dist < std::numeric_limits<double>::epsilon())
                pnt_in_quadtree = true;
        }
        if (!pnt_in_quadtree)
            _pnts.push_back (pnt);
        else
            return false;

        if (_pnts.size () > _max_points_per_leaf)
            splitNode ();
        return true;
    }

    /**
     * This method balances the quadtree, i.e., it will be inserted nodes
     * such that the depth between neighbored leafs is at most one. If you want
     * to create a mesh (for instance with GMSH) you can use this method to
     * improve the mesh quality. The balance method should be used after
     * inserting all points.
     */
    void balance ()
    {
        std::list<QuadTree<POINT>*> leaf_list;
        getLeafs (leaf_list);

        while (!leaf_list.empty())
        {
            QuadTree<POINT>* node (leaf_list.front());
            leaf_list.pop_front ();

            if (node->isLeaf())
                if (needToRefine (node))
                {
                    node->splitNode ();
                    leaf_list.push_back (node->getChild(Quadrant::NE));
                    leaf_list.push_back (node->getChild(Quadrant::NW));
                    leaf_list.push_back (node->getChild(Quadrant::SW));
                    leaf_list.push_back (node->getChild(Quadrant::SE));

                    // check if north neighbor has to be refined
                    QuadTree<POINT>* north_neighbor (node->getNorthNeighbor());
                    if (north_neighbor != nullptr)
                        if (north_neighbor->getDepth() < node->getDepth ())
                            if (north_neighbor->isLeaf())
                                leaf_list.push_back (north_neighbor);

                    // check if west neighbor has to be refined
                    QuadTree<POINT>* west_neighbor (node->getWestNeighbor());
                    if (west_neighbor != nullptr)
                        if (west_neighbor->getDepth() < node->getDepth ())
                            if (west_neighbor->isLeaf())
                                leaf_list.push_back (west_neighbor);

                    // check if south neighbor has to be refined
                    QuadTree<POINT>* south_neighbor (node->getSouthNeighbor());
                    if (south_neighbor != nullptr)
                        if (south_neighbor->getDepth() < node->getDepth ())
                            if (south_neighbor->isLeaf())
                                leaf_list.push_back (south_neighbor);

                    // check if east neighbor has to be refined
                    QuadTree<POINT>* east_neighbor (node->getEastNeighbor());
                    if (east_neighbor != nullptr)
                        if (east_neighbor->getDepth() < node->getDepth ())
                            if (east_neighbor->isLeaf())
                                leaf_list.push_back (east_neighbor);

                }
        }
    }

    /**
     * add all leafs of the quadtree to the list
     * @param leaf_list list of leafs
     */
    void getLeafs (std::list<QuadTree<POINT>*>& leaf_list)
    {
        if (_is_leaf)
            leaf_list.push_back (this);
        else
            for (auto& child : _children)
                child->getLeafs(leaf_list);
    }

    const std::vector<POINT const*>& getPoints () const { return _pnts; }

    void getSquarePoints (POINT& ll, POINT& ur) const
    {
        ll = _ll;
        ur = _ur;
    }

    void getLeaf (const POINT& pnt, POINT& ll, POINT& ur)
    {
        if (this->isLeaf())
        {
            ll = _ll;
            ur = _ur;
        }
        else
        {
            if (pnt[0] <= 0.5 * (_ur[0] + _ll[0])) // WEST
            {
                if (pnt[1] <= 0.5 * (_ur[1] + _ll[1])) // SOUTH
                    _children[static_cast<int>(Quadrant::SW)]->getLeaf (pnt, ll, ur);
                else // NORTH
                    _children[static_cast<int>(Quadrant::NW)]->getLeaf (pnt, ll, ur);
            }
            else // EAST
            {
                if (pnt[1] <= 0.5 * (_ur[1] + _ll[1])) // SOUTH
                    _children[static_cast<int>(Quadrant::SE)]->getLeaf (pnt, ll, ur);
                else // NORTH
                    _children[static_cast<int>(Quadrant::NE)]->getLeaf (pnt, ll, ur);
            }
        }
    }

    QuadTree<POINT> const* getFather ()
    {
        return _father;
    }

    QuadTree<POINT> const* getChild (Quadrant quadrant) const
    {
        return _children[quadrant];
    }

    /**
     * Method calculates the maximum depth of the QuadTree instance. It is used within
     * the method GMSHAdaptiveMeshDensity::getSteinerPoints().
     * @param max_depth (input/output) at the entry max_depth contains the maximum depth up to now
     */
    void getMaxDepth (std::size_t &max_depth) const
    {
        if (max_depth < _depth)
            max_depth = _depth;

        for (auto& child : _children)
        {
            if (child)
            {
                child->getMaxDepth(max_depth);
            }
        }
    }

    /**
     * Method returns the depth of the current QuadTree node.
     * @return the depth of the current QuadTree node
     */
    std::size_t getDepth () const { return _depth; }

private:
    QuadTree<POINT>* getChild (Quadrant quadrant)
    {
        return _children[static_cast<int>(quadrant)];
    }

    bool isLeaf () const { return _is_leaf; }

    bool isChild (QuadTree<POINT> const* const tree, Quadrant quadrant) const
    {
        if (_children[static_cast<int>(quadrant)] == tree) return true;
        return false;
    }

    QuadTree<POINT>* getNorthNeighbor () const
    {
        if (this->_father == nullptr) // root of QuadTree
            return nullptr;

        if (this->_father->isChild (this, Quadrant::SW))
            return this->_father->getChild (Quadrant::NW);
        if (this->_father->isChild (this, Quadrant::SE))
            return this->_father->getChild (Quadrant::NE);

        QuadTree<POINT>* north_neighbor (this->_father->getNorthNeighbor ());
        if (north_neighbor == nullptr)
            return nullptr;
        if (north_neighbor->isLeaf())
            return north_neighbor;

        if (this->_father->isChild (this, Quadrant::NW))
            return north_neighbor->getChild (Quadrant::SW);

        return north_neighbor->getChild(Quadrant::SE);
    }

    QuadTree<POINT>* getSouthNeighbor () const
    {
        if (this->_father == nullptr) // root of QuadTree
            return nullptr;

        if (this->_father->isChild (this, Quadrant::NW))
            return this->_father->getChild (Quadrant::SW);
        if (this->_father->isChild (this, Quadrant::NE))
            return this->_father->getChild (Quadrant::SE);

        QuadTree<POINT>* south_neighbor (this->_father->getSouthNeighbor ());
        if (south_neighbor == nullptr)
            return nullptr;
        if (south_neighbor->isLeaf())
            return south_neighbor;

        if (this->_father->isChild (this, Quadrant::SW))
            return south_neighbor->getChild (Quadrant::NW);

        return south_neighbor->getChild(Quadrant::NE);
    }

    QuadTree<POINT>* getEastNeighbor () const
    {
        if (this->_father == nullptr) // root of QuadTree
            return nullptr;

        if (this->_father->isChild (this, Quadrant::NW))
            return this->_father->getChild (Quadrant::NE);
        if (this->_father->isChild (this, Quadrant::SW))
            return this->_father->getChild (Quadrant::SE);

        QuadTree<POINT>* east_neighbor (this->_father->getEastNeighbor ());
        if (east_neighbor == nullptr)
            return nullptr;
        if (east_neighbor->isLeaf())
            return east_neighbor;

        if (this->_father->isChild (this, Quadrant::SE))
            return east_neighbor->getChild (Quadrant::SW);

        return east_neighbor->getChild(Quadrant::NW);
    }

    QuadTree<POINT>* getWestNeighbor () const
    {
        if (this->_father == nullptr) // root of QuadTree
            return nullptr;

        if (this->_father->isChild (this, Quadrant::NE))
            return this->_father->getChild (Quadrant::NW);
        if (this->_father->isChild (this, Quadrant::SE))
            return this->_father->getChild (Quadrant::SW);

        QuadTree<POINT>* west_neighbor (this->_father->getWestNeighbor ());
        if (west_neighbor == nullptr)
            return nullptr;
        if (west_neighbor->isLeaf())
            return west_neighbor;

        if (this->_father->isChild (this, Quadrant::SW))
            return west_neighbor->getChild (Quadrant::SE);

        return west_neighbor->getChild(Quadrant::NE);
    }

    /**
     * private constructor
     * @param ll lower left point
     * @param ur upper right point
     * @param father father in the tree
     * @param depth depth of the node
     * @param max_points_per_leaf maximum number of points per leaf
     * @return
     */
    QuadTree(POINT ll,
             POINT ur,
             QuadTree* father,
             std::size_t depth,
             std::size_t max_points_per_leaf)
        : _father(father),
          _ll(std::move(ll)),
          _ur(std::move(ur)),
          _depth(depth),
          _is_leaf(true),
          _max_points_per_leaf(max_points_per_leaf)
    {
        // init children
        for (auto& child : _children)
            child = nullptr;
    }

    void splitNode ()
    {
        // create children
        POINT mid_point(_ll);
        mid_point[0] += (_ur[0] - _ll[0]) / 2.0;
        mid_point[1] += (_ur[1] - _ll[1]) / 2.0;
        assert(_children[0] == nullptr);
        _children[0] = new QuadTree<POINT> (mid_point, _ur, this, _depth + 1, _max_points_per_leaf); // north east
        POINT h_ll(mid_point), h_ur(mid_point);
        h_ll[0] = _ll[0];
        h_ur[1] = _ur[1];
        assert(_children[1] == nullptr);
        _children[1] = new QuadTree<POINT> (h_ll, h_ur, this, _depth + 1, _max_points_per_leaf); // north west
        assert(_children[2] == nullptr);
        _children[2] = new QuadTree<POINT> (_ll, mid_point, this, _depth + 1, _max_points_per_leaf); // south west
        h_ll = _ll;
        h_ll[0] = mid_point[0];
        h_ur = _ur;
        h_ur[1] = mid_point[1];
        assert(_children[3] == nullptr);
        _children[3] = new QuadTree<POINT> (h_ll, h_ur, this, _depth + 1, _max_points_per_leaf); // south east

        // distribute points to sub quadtrees
        for (std::size_t j(0); j < _pnts.size(); j++) {
            bool nfound(true);
            for (std::size_t k(0); k < 4 && nfound; k++)
                if (_children[k]->addPoint(_pnts[j])) nfound = false;

        }
        _pnts.clear();
        _is_leaf = false;
    }

    bool needToRefine (QuadTree<POINT>* node)
    {
        QuadTree<POINT>* north_neighbor (node->getNorthNeighbor ());
        if (north_neighbor != nullptr)
        {
            if (north_neighbor->getDepth() == node->getDepth())
                if (!north_neighbor->isLeaf ())
                {
                    if (!(north_neighbor->getChild(Quadrant::SW))->isLeaf())
                        return true;
                    if (!(north_neighbor->getChild(Quadrant::SE))->isLeaf())
                        return true;
                }
        }

        QuadTree<POINT>* west_neighbor (node->getWestNeighbor ());
        if (west_neighbor != nullptr)
        {
            if (west_neighbor->getDepth() == node->getDepth())
                if (!west_neighbor->isLeaf ())
                {
                    if (!(west_neighbor->getChild(Quadrant::SE))->isLeaf())
                        return true;
                    if (!(west_neighbor->getChild(Quadrant::NE))->isLeaf())
                        return true;
                }
        }

        QuadTree<POINT>* south_neighbor (node->getSouthNeighbor ());
        if (south_neighbor != nullptr)
        {
            if (south_neighbor->getDepth() == node->getDepth())
                if (!south_neighbor->isLeaf())
                {
                    if (!(south_neighbor->getChild(Quadrant::NE))->isLeaf())
                        return true;
                    if (!(south_neighbor->getChild(Quadrant::NW))->isLeaf())
                        return true;
                }
        }

        QuadTree<POINT>* east_neighbor (node->getEastNeighbor ());
        if (east_neighbor != nullptr)
        {
            if (east_neighbor->getDepth() == node->getDepth())
                if (!east_neighbor->isLeaf ())
                {
                    if (!(east_neighbor->getChild(Quadrant::NW))->isLeaf())
                        return true;
                    if (!(east_neighbor->getChild(Quadrant::SW))->isLeaf())
                        return true;
                }
        }
        return false;
    }

    QuadTree<POINT>* _father;
    /**
     * children are sorted:
     *   _children[0] is north east child
     *   _children[1] is north west child
     *   _children[2] is south west child
     *   _children[3] is south east child
     */
    QuadTree<POINT>* _children[4];
    /**
     * lower left point of the square
     */
    POINT _ll;
    /**
     * upper right point of the square
     */
    POINT _ur;
    std::size_t _depth;
    std::vector<POINT const*> _pnts;
    bool _is_leaf;
    /**
     * maximum number of points per leaf
     */
    const std::size_t _max_points_per_leaf;
};
}
