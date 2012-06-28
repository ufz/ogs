/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file QuadTree.h
 *
 * Created on 2010-11-09 by Thomas Fischer
 */

#ifndef QUADTREE_H_
#define QUADTREE_H_

namespace GeoLib {

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
template <typename POINT> class QuadTree {
public:
	enum Quadrant {
		NE = 0, //!< north east
		NW, //!< north west
		SW, //!< south west
		SE  //!< south east
	};
	/**
	 * This is the constructor for class QuadTree. It takes two points
	 * (lower left and the upper right points).
	 * @param ll lower left point of the square
	 * @param ur upper right point of the square
	 * @param max_points_per_node maximum number of points per node, if the
	 * maximum number of points per node is reached and one point more should
	 * be added the node will be split and the points are distributed to the
	 * childs of the node
	 */
	QuadTree(POINT const& ll, POINT const& ur, size_t max_points_per_node) :
		_father (NULL), _ll (ll), _ur (ur), _depth (0), _is_leaf (true),
		_max_points_per_node (max_points_per_node)
	{
		assert (_max_points_per_node > 0);

		// init childs
		for (size_t k(0); k<4; k++) {
			_childs[k] = NULL;
		}

		if ((_ur[0] - _ll[0]) > (_ur[1] - _ll[1])) {
			_ur[1] = _ll[1] + _ur[0] - _ll[0];
		} else {
			_ur[0] = _ll[0] + _ur[1] - _ll[1];
		}
//#ifndef NDEBUG
//		std::cerr << "lower left: " << _ll << ", upper right: " << _ur << ", depth " << _depth << std::endl;
//#endif
	}

	/**
	 * destructor
	 */
	~QuadTree()
	{
		if (_is_leaf) {
			for (size_t k(0); k<4; k++) {
				delete _childs[k];
			}
		}
	}

	/**
	 * This method adds the given point to the quadtree. If necessary,
	 * the quadtree will be extended.
	 * @param pnt the point
	 * @return If the point can be inserted the method returns true, else false.
	 */
	bool addPoint (POINT * pnt)
	{
		if ((*pnt)[0] < _ll[0]) return false;
		if ((*pnt)[0] > _ur[0]) return false;
		if ((*pnt)[1] < _ll[1]) return false;
		if ((*pnt)[1] > _ur[1]) return false;

		if (!_is_leaf) {
			for (size_t k(0); k<4; k++) {
				if (_childs[k]->addPoint (pnt))
					return true;
			}
		}

		// check if point is already in quadtree
		bool pnt_in_quadtree (false);
		double equal_pnt_dist (MathLib::fastpow(2, _depth) * fabs(_ll[0] - _ur[0]) * 1e-6);
		for (size_t k(0); k<_pnts.size() && !pnt_in_quadtree; k++) {
			const double sqr_dist (MathLib::sqrDist(_pnts[k]->getCoords(), pnt->getCoords()));
			if (sqr_dist < equal_pnt_dist) {
				pnt_in_quadtree = true;
			}
		}
		if (!pnt_in_quadtree) {
			_pnts.push_back (pnt);
		} else {
			return false;
		}

		if (_pnts.size () > _max_points_per_node) {
			splitNode ();
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

		while (!leaf_list.empty()) {
			QuadTree<POINT>* node (leaf_list.front());
			leaf_list.pop_front ();

			if (node->isLeaf()) {
				if (needToRefine (node)) {
					node->splitNode ();
					leaf_list.push_back (node->getChild(NE));
					leaf_list.push_back (node->getChild(NW));
					leaf_list.push_back (node->getChild(SW));
					leaf_list.push_back (node->getChild(SE));

					// check if north neighbor has to be refined
					QuadTree<POINT>* north_neighbor (node->getNorthNeighbor());
					if (north_neighbor != NULL) {
						if (north_neighbor->getDepth() < node->getDepth ()) {
							if (north_neighbor->isLeaf()) {
								leaf_list.push_back (north_neighbor);
							}
						}
					}

					// check if west neighbor has to be refined
					QuadTree<POINT>* west_neighbor (node->getWestNeighbor());
					if (west_neighbor != NULL) {
						if (west_neighbor->getDepth() < node->getDepth ()) {
							if (west_neighbor->isLeaf()) {
								leaf_list.push_back (west_neighbor);
							}
						}
					}

					// check if south neighbor has to be refined
					QuadTree<POINT>* south_neighbor (node->getSouthNeighbor());
					if (south_neighbor != NULL) {
						if (south_neighbor->getDepth() < node->getDepth ()) {
							if (south_neighbor->isLeaf()) {
								leaf_list.push_back (south_neighbor);
							}
						}
					}

					// check if east neighbor has to be refined
					QuadTree<POINT>* east_neighbor (node->getEastNeighbor());
					if (east_neighbor != NULL) {
						if (east_neighbor->getDepth() < node->getDepth ()) {
							if (east_neighbor->isLeaf()) {
								leaf_list.push_back (east_neighbor);
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
		if (_is_leaf) {
			leaf_list.push_back (this);
		} else {
			for (size_t k(0); k<4; k++) {
				_childs[k]->getLeafs (leaf_list);
			}
		}
	}

	const std::vector<POINT*>& getPoints () const { return _pnts; }

	void getSquarePoints (POINT& ll, POINT& ur) const
	{
		ll = _ll;
		ur = _ur;
	}

	void getLeaf (const POINT& pnt, POINT& ll, POINT& ur)
	{
		if (this->isLeaf()) {
			ll = _ll;
			ur = _ur;
		} else {
			if (pnt[0] <= 0.5*(_ur[0]+_ll[0])) { // WEST
				if (pnt[1] <= 0.5*(_ur[1]+_ll[1])) { // SOUTH
					_childs[SW]->getLeaf (pnt, ll, ur);
				} else { // NORTH
					_childs[NW]->getLeaf (pnt, ll, ur);
				}
			} else { // EAST
				if (pnt[1] <= 0.5*(_ur[1]+_ll[1])) { // SOUTH
					_childs[SE]->getLeaf (pnt, ll, ur);
				} else { // NORTH
					_childs[NE]->getLeaf (pnt, ll, ur);
				}
			}
		}
	}

	void getQuadTree (std::vector<POINT*>& pnts, std::vector<GeoLib::Polyline*>& plys) const
	{
		size_t pnt_pos (pnts.size());
		pnts.push_back (new POINT (_ll));
		pnts.push_back (new POINT (_ur[0], _ll[1], _ll[2]));
		pnts.push_back (new POINT (_ur));
		pnts.push_back (new POINT (_ll[0], _ur[1], _ll[2]));

		if (_father == NULL) {
			size_t ply_pos (plys.size());
			plys.push_back (new Polyline (pnts));
			for (size_t i(0); i<4; i++)
				plys[ply_pos]->addPoint (pnt_pos+i);
			plys[ply_pos]->addPoint (pnt_pos);
		}

		if (! _is_leaf) {
			for (size_t i(0); i<4; i++) {
				_childs[i]->getQuadTree (pnts, plys);
			}
		}
	}

	QuadTree<POINT> const * getFather ()
	{
		return _father;
	}

	QuadTree<POINT> const * getChild (Quadrant quadrant) const
	{
		return _childs[quadrant];
	}


private:
	QuadTree<POINT> * getChild (Quadrant quadrant)
	{
		return _childs[quadrant];
	}

	bool isLeaf () const { return _is_leaf; }

	bool isChild (QuadTree<POINT> const * const tree, Quadrant quadrant) const
	{
		if (_childs[quadrant] == tree) return true;
		return false;
	}

	QuadTree<POINT>* getNorthNeighbor () const
	{
		if (this->_father == NULL) { // root of QuadTree
			return NULL;
		}

		if (this->_father->isChild (this, SW))
			return this->_father->getChild (NW);
		if (this->_father->isChild (this, SE))
			return this->_father->getChild (NE);

		QuadTree<POINT>* north_neighbor (this->_father->getNorthNeighbor ());
		if (north_neighbor == NULL)
			return NULL;
		if (north_neighbor->isLeaf())
			return north_neighbor;

		if (this->_father->isChild (this, NW))
			return north_neighbor->getChild (SW);
		else
			return north_neighbor->getChild (SE);
	}

	QuadTree<POINT>* getSouthNeighbor () const
	{
		if (this->_father == NULL) { // root of QuadTree
			return NULL;
		}

		if (this->_father->isChild (this, NW))
			return this->_father->getChild (SW);
		if (this->_father->isChild (this, NE))
			return this->_father->getChild (SE);

		QuadTree<POINT>* south_neighbor (this->_father->getSouthNeighbor ());
		if (south_neighbor == NULL)
			return NULL;
		if (south_neighbor->isLeaf())
			return south_neighbor;

		if (this->_father->isChild (this, SW))
			return south_neighbor->getChild (NW);
		else
			return south_neighbor->getChild (NE);
	}

	QuadTree<POINT>* getEastNeighbor () const
	{
		if (this->_father == NULL) { // root of QuadTree
			return NULL;
		}

		if (this->_father->isChild (this, NW))
			return this->_father->getChild (NE);
		if (this->_father->isChild (this, SW))
			return this->_father->getChild (SE);

		QuadTree<POINT>* east_neighbor (this->_father->getEastNeighbor ());
		if (east_neighbor == NULL)
			return NULL;
		if (east_neighbor->isLeaf())
			return east_neighbor;

		if (this->_father->isChild (this, SE))
			return east_neighbor->getChild (SW);
		else
			return east_neighbor->getChild (NW);
	}

	QuadTree<POINT>* getWestNeighbor () const
	{
		if (this->_father == NULL) { // root of QuadTree
			return NULL;
		}

		if (this->_father->isChild (this, NE))
			return this->_father->getChild (NW);
		if (this->_father->isChild (this, SE))
			return this->_father->getChild (SW);

		QuadTree<POINT>* west_neighbor (this->_father->getWestNeighbor ());
		if (west_neighbor == NULL)
			return NULL;
		if (west_neighbor->isLeaf())
			return west_neighbor;

		if (this->_father->isChild (this, SW))
			return west_neighbor->getChild (SE);
		else
			return west_neighbor->getChild (NE);
	}

	size_t getDepth () const { return _depth; }

	/**
	 * private constructor
	 * @param ll lower left point
	 * @param ur upper right point
	 * @param father father in the tree
	 * @param depth depth of the node
	 * @param max_points_per_node maximum number of points per node, if the
	 * maximum number of points per node is reached and one point more should
	 * be added the node will be split and the points are distributed to the
	 * @return
	 */
	QuadTree (POINT const& ll, POINT const& ur, QuadTree* father, size_t depth, size_t max_points_per_node) :
		_father (father), _ll (ll), _ur (ur), _depth (depth), _is_leaf (true),
		_max_points_per_node (max_points_per_node)
	{
		// init childs
		for (size_t k(0); k<4; k++) {
			_childs[k] = NULL;
		}

//#ifndef NDEBUG
//		std::cerr << "lower left: " << _ll << ", upper right: " << _ur << ", depth: " << _depth << std::endl;
//#endif
	}

	void splitNode ()
	{
		// create childs
		POINT mid_point(_ll);
		mid_point[0] += (_ur[0] - _ll[0]) / 2.0;
		mid_point[1] += (_ur[1] - _ll[1]) / 2.0;
		_childs[0] = new QuadTree<POINT> (mid_point, _ur, this, _depth + 1, _max_points_per_node); // north east
		POINT h_ll(mid_point), h_ur(mid_point);
		h_ll[0] = _ll[0];
		h_ur[1] = _ur[1];
		_childs[1] = new QuadTree<POINT> (h_ll, h_ur, this, _depth + 1, _max_points_per_node); // north west
		_childs[2] = new QuadTree<POINT> (_ll, mid_point, this, _depth + 1, _max_points_per_node); // south west
		h_ll = _ll;
		h_ll[0] = mid_point[0];
		h_ur = _ur;
		h_ur[1] = mid_point[1];
		_childs[3] = new QuadTree<POINT> (h_ll, h_ur, 	this, _depth + 1, _max_points_per_node); // south east

		// distribute points to sub quadtrees
		for (size_t j(0); j < _pnts.size(); j++) {
			bool nfound(true);
			for (size_t k(0); k < 4 && nfound; k++) {
				if (_childs[k]->addPoint(_pnts[j])) {
					nfound = false;
				}
			}
		}
		_pnts.clear();
		_is_leaf = false;
	}

	bool needToRefine (QuadTree<POINT>* node)
	{
		QuadTree<POINT>* north_neighbor (node->getNorthNeighbor ());
		if (north_neighbor != NULL) {
			if (north_neighbor->getDepth() == node->getDepth()) {
				if (! north_neighbor->isLeaf ()) {
					if (! (north_neighbor->getChild(SW))->isLeaf()) {
						return true;
					}
					if (! (north_neighbor->getChild(SE))->isLeaf()) {
						return true;
					}
				}
			}
		}

		QuadTree<POINT>* west_neighbor (node->getWestNeighbor ());
		if (west_neighbor != NULL) {
			if (west_neighbor->getDepth() == node->getDepth()) {
				if (! west_neighbor->isLeaf ()) {
					if (! (west_neighbor->getChild(SE))->isLeaf()) {
						return true;
					}
					if (! (west_neighbor->getChild(NE))->isLeaf()) {
						return true;
					}
				}
			}
		}

		QuadTree<POINT>* south_neighbor (node->getSouthNeighbor ());
		if (south_neighbor != NULL) {
			if (south_neighbor->getDepth() == node->getDepth()) {
				if (!south_neighbor->isLeaf()) {
					if (!(south_neighbor->getChild(NE))->isLeaf()) {
						return true;
					}
					if (!(south_neighbor->getChild(NW))->isLeaf()) {
						return true;
					}
				}
			}
		}

		QuadTree<POINT>* east_neighbor (node->getEastNeighbor ());
		if (east_neighbor != NULL) {
			if (east_neighbor->getDepth() == node->getDepth()) {
				if (! east_neighbor->isLeaf ()) {
					if (! (east_neighbor->getChild(NW))->isLeaf()) {
						return true;
					}
					if (! (east_neighbor->getChild(SW))->isLeaf()) {
						return true;
					}
				}
			}
		}
		return false;
	}

	QuadTree<POINT>* _father;
	/**
	 * childs are sorted:
	 *   _childs[0] is north east child
	 *   _childs[1] is north west child
	 *   _childs[2] is south west child
	 *   _childs[3] is south east child
	 */
	QuadTree<POINT>* _childs[4];
	/**
	 * lower left point of the square
	 */
	POINT _ll;
	/**
	 * upper right point of the square
	 */
	POINT _ur;
	size_t _depth;
	std::vector<POINT *> _pnts;
	bool _is_leaf;
	/**
	 * maximum number of points per leaf
	 */
	const size_t _max_points_per_node;
};

}

#endif /* QUADTREE_H_ */
