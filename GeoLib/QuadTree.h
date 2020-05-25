/**
 * \file
 * \author Thomas Fischer
 * \date   2010-11-09
 * \brief  Definition of the QuadTree class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <limits>
#include <list>
#include <utility>

#include "BaseLib/Logging.h"

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
        : father_(nullptr),
          ll_(std::move(ll)),
          ur_(std::move(ur)),
          depth_(0),
          is_leaf_(true),
          max_points_per_leaf_(max_points_per_leaf)
    {
        assert (max_points_per_leaf_ > 0);

        // init children
        for (auto& child : children_)
        {
            child = nullptr;
        }

        if ((ur_[0] - ll_[0]) > (ur_[1] - ll_[1]))
        {
            ur_[1] = ll_[1] + ur_[0] - ll_[0];
        }
        else
        {
            ur_[0] = ll_[0] + ur_[1] - ll_[1];
        }

        DBUG(
            "QuadTree(): lower left: ({:f},{:f},{:f}), upper right: "
            "({:f},{:f},{:f}), depth: {:d}",
            ll_[0],
            ll_[1],
            ll_[2],
            ur_[0],
            ur_[1],
            ur_[2],
            depth_);
    }

    /**
     * destructor
     */
    ~QuadTree()
    {
        for (auto& child : children_)
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
        if ((*pnt)[0] < ll_[0])
        {
            return false;
        }
        if ((*pnt)[0] >= ur_[0])
        {
            return false;
        }
        if ((*pnt)[1] < ll_[1])
        {
            return false;
        }
        if ((*pnt)[1] >= ur_[1])
        {
            return false;
        }

        if (!is_leaf_) {
            for (auto& child : children_)
            {
                if (child->addPoint(pnt))
                {
                    return true;
                }
            }
            return false;
        }

        // check if point is already in quadtree
        bool pnt_in_quadtree (false);
        for (std::size_t k(0); k < pnts_.size() && !pnt_in_quadtree; k++) {
            const double v0((*(pnts_[k]))[0] - (*pnt)[0]);
            const double v1((*(pnts_[k]))[1] - (*pnt)[1]);
            const double sqr_dist (v0*v0 + v1*v1);
            if (sqr_dist < std::numeric_limits<double>::epsilon())
            {
                pnt_in_quadtree = true;
            }
        }
        if (!pnt_in_quadtree)
        {
            pnts_.push_back (pnt);
        }
        else
        {
            return false;
        }

        if (pnts_.size() > max_points_per_leaf_)
        {
            splitNode();
        }
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
            {
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
                    {
                        if (north_neighbor->getDepth() < node->getDepth())
                        {
                            if (north_neighbor->isLeaf())
                            {
                                leaf_list.push_back(north_neighbor);
                            }
                        }
                    }

                    // check if west neighbor has to be refined
                    QuadTree<POINT>* west_neighbor (node->getWestNeighbor());
                    if (west_neighbor != nullptr)
                    {
                        if (west_neighbor->getDepth() < node->getDepth())
                        {
                            if (west_neighbor->isLeaf())
                            {
                                leaf_list.push_back(west_neighbor);
                            }
                        }
                    }

                    // check if south neighbor has to be refined
                    QuadTree<POINT>* south_neighbor (node->getSouthNeighbor());
                    if (south_neighbor != nullptr)
                    {
                        if (south_neighbor->getDepth() < node->getDepth())
                        {
                            if (south_neighbor->isLeaf())
                            {
                                leaf_list.push_back(south_neighbor);
                            }
                        }
                    }

                    // check if east neighbor has to be refined
                    QuadTree<POINT>* east_neighbor (node->getEastNeighbor());
                    if (east_neighbor != nullptr)
                    {
                        if (east_neighbor->getDepth() < node->getDepth())
                        {
                            if (east_neighbor->isLeaf())
                            {
                                leaf_list.push_back(east_neighbor);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * add all leafs of the quadtree to the list
     * @param leaf_list list of leafs
     */
    void getLeafs (std::list<QuadTree<POINT>*>& leaf_list)
    {
        if (is_leaf_)
        {
            leaf_list.push_back (this);
        }
        else
        {
            for (auto& child : children_)
            {
                child->getLeafs(leaf_list);
            }
        }
    }

    const std::vector<POINT const*>& getPoints () const { return pnts_; }

    void getSquarePoints (POINT& ll, POINT& ur) const
    {
        ll = ll_;
        ur = ur_;
    }

    void getLeaf (const POINT& pnt, POINT& ll, POINT& ur)
    {
        if (this->isLeaf())
        {
            ll = ll_;
            ur = ur_;
        }
        else
        {
            if (pnt[0] <= 0.5 * (ur_[0] + ll_[0])) // WEST
            {
                if (pnt[1] <= 0.5 * (ur_[1] + ll_[1]))
                {  // SOUTH
                    children_[static_cast<int>(Quadrant::SW)]->getLeaf (pnt, ll, ur);
                }
                else
                {  // NORTH
                    children_[static_cast<int>(Quadrant::NW)]->getLeaf(
                        pnt, ll, ur);
                }
            }
            else // EAST
            {
                if (pnt[1] <= 0.5 * (ur_[1] + ll_[1]))
                {  // SOUTH
                    children_[static_cast<int>(Quadrant::SE)]->getLeaf (pnt, ll, ur);
                }
                else
                {  // NORTH
                    children_[static_cast<int>(Quadrant::NE)]->getLeaf(
                        pnt, ll, ur);
                }
            }
        }
    }

    QuadTree<POINT> const* getFather ()
    {
        return father_;
    }

    QuadTree<POINT> const* getChild (Quadrant quadrant) const
    {
        return children_[quadrant];
    }

    /**
     * Method calculates the maximum depth of the QuadTree instance. It is used within
     * the method GMSHAdaptiveMeshDensity::getSteinerPoints().
     * @param max_depth (input/output) at the entry max_depth contains the maximum depth up to now
     */
    void getMaxDepth (std::size_t &max_depth) const
    {
        if (max_depth < depth_)
        {
            max_depth = depth_;
        }

        for (auto& child : children_)
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
    std::size_t getDepth () const { return depth_; }

private:
    QuadTree<POINT>* getChild (Quadrant quadrant)
    {
        return children_[static_cast<int>(quadrant)];
    }

    bool isLeaf () const { return is_leaf_; }

    bool isChild (QuadTree<POINT> const* const tree, Quadrant quadrant) const
    {
        return children_[static_cast<int>(quadrant)] == tree;
    }

    QuadTree<POINT>* getNorthNeighbor () const
    {
        if (this->father_ == nullptr)
        {  // root of QuadTree
            return nullptr;
        }

        if (this->father_->isChild(this, Quadrant::SW))
        {
            return this->father_->getChild(Quadrant::NW);
        }
        if (this->father_->isChild(this, Quadrant::SE))
        {
            return this->father_->getChild(Quadrant::NE);
        }

        QuadTree<POINT>* north_neighbor (this->father_->getNorthNeighbor ());
        if (north_neighbor == nullptr)
        {
            return nullptr;
        }
        if (north_neighbor->isLeaf())
        {
            return north_neighbor;
        }

        if (this->father_->isChild(this, Quadrant::NW))
        {
            return north_neighbor->getChild(Quadrant::SW);
        }

        return north_neighbor->getChild(Quadrant::SE);
    }

    QuadTree<POINT>* getSouthNeighbor () const
    {
        if (this->father_ == nullptr)
        {  // root of QuadTree
            return nullptr;
        }

        if (this->father_->isChild(this, Quadrant::NW))
        {
            return this->father_->getChild(Quadrant::SW);
        }
        if (this->father_->isChild(this, Quadrant::NE))
        {
            return this->father_->getChild(Quadrant::SE);
        }

        QuadTree<POINT>* south_neighbor (this->father_->getSouthNeighbor ());
        if (south_neighbor == nullptr)
        {
            return nullptr;
        }
        if (south_neighbor->isLeaf())
        {
            return south_neighbor;
        }

        if (this->father_->isChild(this, Quadrant::SW))
        {
            return south_neighbor->getChild(Quadrant::NW);
        }

        return south_neighbor->getChild(Quadrant::NE);
    }

    QuadTree<POINT>* getEastNeighbor () const
    {
        if (this->father_ == nullptr)
        {  // root of QuadTree
            return nullptr;
        }

        if (this->father_->isChild(this, Quadrant::NW))
        {
            return this->father_->getChild(Quadrant::NE);
        }
        if (this->father_->isChild(this, Quadrant::SW))
        {
            return this->father_->getChild(Quadrant::SE);
        }

        QuadTree<POINT>* east_neighbor (this->father_->getEastNeighbor ());
        if (east_neighbor == nullptr)
        {
            return nullptr;
        }
        if (east_neighbor->isLeaf())
        {
            return east_neighbor;
        }

        if (this->father_->isChild(this, Quadrant::SE))
        {
            return east_neighbor->getChild(Quadrant::SW);
        }

        return east_neighbor->getChild(Quadrant::NW);
    }

    QuadTree<POINT>* getWestNeighbor () const
    {
        if (this->father_ == nullptr)
        {  // root of QuadTree
            return nullptr;
        }

        if (this->father_->isChild(this, Quadrant::NE))
        {
            return this->father_->getChild(Quadrant::NW);
        }
        if (this->father_->isChild(this, Quadrant::SE))
        {
            return this->father_->getChild(Quadrant::SW);
        }

        QuadTree<POINT>* west_neighbor (this->father_->getWestNeighbor ());
        if (west_neighbor == nullptr)
        {
            return nullptr;
        }
        if (west_neighbor->isLeaf())
        {
            return west_neighbor;
        }

        if (this->father_->isChild(this, Quadrant::SW))
        {
            return west_neighbor->getChild(Quadrant::SE);
        }

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
        : father_(father),
          ll_(std::move(ll)),
          ur_(std::move(ur)),
          depth_(depth),
          is_leaf_(true),
          max_points_per_leaf_(max_points_per_leaf)
    {
        // init children
        for (auto& child : children_)
        {
            child = nullptr;
        }
    }

    void splitNode ()
    {
        // create children
        POINT mid_point(ll_);
        mid_point[0] += (ur_[0] - ll_[0]) / 2.0;
        mid_point[1] += (ur_[1] - ll_[1]) / 2.0;
        assert(children_[0] == nullptr);
        children_[0] = new QuadTree<POINT> (mid_point, ur_, this, depth_ + 1, max_points_per_leaf_); // north east
        POINT h_ll(mid_point), h_ur(mid_point);
        h_ll[0] = ll_[0];
        h_ur[1] = ur_[1];
        assert(children_[1] == nullptr);
        children_[1] = new QuadTree<POINT> (h_ll, h_ur, this, depth_ + 1, max_points_per_leaf_); // north west
        assert(children_[2] == nullptr);
        children_[2] = new QuadTree<POINT> (ll_, mid_point, this, depth_ + 1, max_points_per_leaf_); // south west
        h_ll = ll_;
        h_ll[0] = mid_point[0];
        h_ur = ur_;
        h_ur[1] = mid_point[1];
        assert(children_[3] == nullptr);
        children_[3] = new QuadTree<POINT> (h_ll, h_ur, this, depth_ + 1, max_points_per_leaf_); // south east

        // distribute points to sub quadtrees
        for (std::size_t j(0); j < pnts_.size(); j++) {
            bool nfound(true);
            for (std::size_t k(0); k < 4 && nfound; k++)
            {
                if (children_[k]->addPoint(pnts_[j]))
                {
                    nfound = false;
                }
            }
        }
        pnts_.clear();
        is_leaf_ = false;
    }

    bool needToRefine (QuadTree<POINT>* node)
    {
        QuadTree<POINT>* north_neighbor (node->getNorthNeighbor ());
        if (north_neighbor != nullptr)
        {
            if (north_neighbor->getDepth() == node->getDepth())
            {
                if (!north_neighbor->isLeaf ())
                {
                    if (!(north_neighbor->getChild(Quadrant::SW))->isLeaf())
                    {
                        return true;
                    }
                    if (!(north_neighbor->getChild(Quadrant::SE))->isLeaf())
                    {
                        return true;
                    }
                }
            }
        }

        QuadTree<POINT>* west_neighbor (node->getWestNeighbor ());
        if (west_neighbor != nullptr)
        {
            if (west_neighbor->getDepth() == node->getDepth())
            {
                if (!west_neighbor->isLeaf ())
                {
                    if (!(west_neighbor->getChild(Quadrant::SE))->isLeaf())
                    {
                        return true;
                    }
                    if (!(west_neighbor->getChild(Quadrant::NE))->isLeaf())
                    {
                        return true;
                    }
                }
            }
        }

        QuadTree<POINT>* south_neighbor (node->getSouthNeighbor ());
        if (south_neighbor != nullptr)
        {
            if (south_neighbor->getDepth() == node->getDepth())
            {
                if (!south_neighbor->isLeaf())
                {
                    if (!(south_neighbor->getChild(Quadrant::NE))->isLeaf())
                    {
                        return true;
                    }
                    if (!(south_neighbor->getChild(Quadrant::NW))->isLeaf())
                    {
                        return true;
                    }
                }
            }
        }

        QuadTree<POINT>* east_neighbor (node->getEastNeighbor ());
        if (east_neighbor != nullptr)
        {
            if (east_neighbor->getDepth() == node->getDepth())
            {
                if (!east_neighbor->isLeaf ())
                {
                    if (!(east_neighbor->getChild(Quadrant::NW))->isLeaf())
                    {
                        return true;
                    }
                    if (!(east_neighbor->getChild(Quadrant::SW))->isLeaf())
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    QuadTree<POINT>* father_;
    /**
     * children are sorted:
     *   children_[0] is north east child
     *   children_[1] is north west child
     *   children_[2] is south west child
     *   children_[3] is south east child
     */
    QuadTree<POINT>* children_[4];
    /**
     * lower left point of the square
     */
    POINT ll_;
    /**
     * upper right point of the square
     */
    POINT ur_;
    std::size_t depth_;
    std::vector<POINT const*> pnts_;
    bool is_leaf_;
    /**
     * maximum number of points per leaf
     */
    const std::size_t max_points_per_leaf_;
};
}  // namespace GeoLib
