/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file  OctTree.h
 *
 * Created on 2012-02-27 by Thomas Fischer
 */

#ifndef OCTTREE_H_
#define OCTTREE_H_

namespace GeoLib {

template <typename POINT> class OctTree {
public:

	static OctTree<POINT>* createOctTree(POINT & ll, POINT & ur, std::size_t max_points_per_node)
	{
		const double dx(ur[0] - ll[0]);
		const double dy(ur[1] - ll[1]);
		const double dz(ur[2] - ll[2]);

		if (dx >= dy && dx >= dz) {
			ll[1] -= (dx-dy)/2.0;
			ur[1] += (dx-dy)/2.0;
			ll[2] -= (dx-dz)/2.0;
			ur[2] += (dx-dz)/2.0;
		} else {
			if (dy >= dx && dy >= dz) {
				ll[0] -= (dy-dx)/2.0;
				ur[0] += (dy-dx)/2.0;
				ll[2] -= (dy-dz)/2.0;
				ur[2] += (dy-dz)/2.0;
			} else {
				ll[0] -= (dz-dx)/2.0;
				ur[0] += (dz-dx)/2.0;
				ll[1] -= (dz-dy)/2.0;
				ur[1] += (dz-dy)/2.0;
			}
		}

		OctTree<POINT>::_max_points_per_node = max_points_per_node;
		return new OctTree<POINT>(ll, ur);
	}

	virtual ~OctTree()
	{
		for (std::size_t k(0); k < 8; k++)
			delete _childs[k];
	}

	/**
	 * This method adds the given point to the OctTree. If necessary,
	 * the OctTree will be extended.
	 * @param pnt the point
	 * @return If the point can be inserted the method returns true, else false.
	 */
	bool addPoint (POINT* pnt)
	{
		if ((*pnt)[0] < _ll[0]) return false;
		if ((*pnt)[0] > _ur[0]) return false;
		if ((*pnt)[1] < _ll[1]) return false;
		if ((*pnt)[1] > _ur[1]) return false;
		if ((*pnt)[2] < _ll[2]) return false;
		if ((*pnt)[2] > _ur[2]) return false;

		if (!_is_leaf) {
			for (std::size_t k(0); k < 8; k++) {
				if (_childs[k]->addPoint (pnt)) {
					return true;
				}
			}
		}

		// check if point is already in OctTree
		bool pnt_in_tree (false);
		for (std::size_t k(0); k < _pnts.size() && !pnt_in_tree; k++) {
			const double sqr_dist (MathLib::sqrDist( (_pnts[k])->getCoords(), pnt->getCoords() ));
			if (sqr_dist < std::numeric_limits<double>::epsilon())
				pnt_in_tree = true;
		}
		if (!pnt_in_tree)
			_pnts.push_back (pnt);
		else
			return false;

		if (_pnts.size () > OctTree<POINT>::_max_points_per_node)
			splitNode ();
		return true;
	}

	/**
	 * range query - returns all points inside the range (min[0], max[0]) x (min[1], max[1]) x (min[2], max[2])
	 * @param min
	 * @param max
	 * @param pnts
	 */
	void getPointsInRange(POINT const& min, POINT const& max, std::vector<POINT*> &pnts) const
	{
		if (_ur[0] < min[0]) return;
		if (_ur[1] < min[1]) return;
		if (_ur[2] < min[2]) return;

		if (max[0] < _ll[0]) return;
		if (max[1] < _ll[1]) return;
		if (max[2] < _ll[2]) return;

		if (_is_leaf) {
			typename std::vector<POINT*>::const_iterator it;
			for (it = (_pnts.begin()); it != _pnts.end(); it++) {
				pnts.push_back(*it);
			}
		} else {
			for (std::size_t k(0); k<8; k++) {
				_childs[k]->getPointsInRange(min, max, pnts);
			}
		}
	}

private:
	enum OctTreeQuadrant {
		NEL = 0, //!< north east lower
		NWL, //!< north west lower
		SWL, //!< south west lower
		SEL, //!< south east lower
		NEU, //!< south west upper
		NWU, //!< south west upper
		SWU, //!< south west upper
		SEU //!< south east upper
	};

	/**
	 * private constructor
	 * @param ll lower left point
	 * @param ur upper right point
	 * @return
	 */
	OctTree (POINT const& ll, POINT const& ur) :
		_ll (ll), _ur (ur), _is_leaf (true)
	{
		// init childs
		for (std::size_t k(0); k < 8; k++)
			_childs[k] = NULL;
	}

	void splitNode ()
	{
		const double x_mid((_ur[0] + _ll[0]) / 2.0);
		const double y_mid((_ur[1] + _ll[1]) / 2.0);
		const double z_mid((_ur[2] + _ll[2]) / 2.0);
		POINT p0(x_mid, y_mid, _ll[2]), p1(_ur[0], _ur[1], z_mid);

		// create child NEL
		_childs[NEL] = new OctTree<POINT> (p0, p1);

		// create child NWL
		p0[0] = _ll[0];
		p1[0] = x_mid;
		_childs[NWL] = new OctTree<POINT> (p0, p1);

		// create child SWL
		p0[1] = _ll[1];
		p1[1] = y_mid;
		_childs[SWL] = new OctTree<POINT> (_ll, p1);

		// create child NEU
		_childs[NEU] = new OctTree<POINT> (p1, _ur);

		// create child SEL
		p0[0] = x_mid;
		p1[0] = _ur[0];
		_childs[SEL] = new OctTree<POINT> (p0, p1);

		// create child NWU
		p0[0] = _ll[0];
		p0[1] = y_mid;
		p0[2] = z_mid;
		p1[0] = x_mid;
		p1[1] = _ur[1];
		p1[2] = _ur[2];
		_childs[NWU] = new OctTree<POINT> (p0, p1);

		// create child SWU
		p0[1] = _ll[1];
		p1[1] = y_mid;
		_childs[SWU] = new OctTree<POINT> (p0, p1);

		// create child SEU
		p0[0] = x_mid;
		p1[0] = _ur[0];
		p1[1] = y_mid;
		p1[2] = _ur[2];
		_childs[SEU] = new OctTree<POINT> (p0, p1);

		// distribute points to sub quadtrees
		const std::size_t n_pnts(_pnts.size());
		for (std::size_t j(0); j < n_pnts; j++) {
			bool nfound(true);
			for (std::size_t k(0); k < 8 && nfound; k++) {
				if (_childs[k]->addPoint(_pnts[j])) {
					nfound = false;
				}
			}
		}
		_pnts.clear();
		_is_leaf = false;
	}

	/**
	 * childs are sorted:
	 *   _childs[0] is north east lower child
	 *   _childs[1] is north west lower child
	 *   _childs[2] is south west lower child
	 *   _childs[3] is south east lower child
	 *   _childs[4] is north east upper child
	 *   _childs[5] is north west upper child
	 *   _childs[6] is south west upper child
	 *   _childs[7] is south east upper child
	 */
	OctTree<POINT>* _childs[8];
	/**
	 * lower left front face point of the cube
	 */
	POINT const _ll;
	/**
	 * upper right back face point of the cube
	 */
	POINT const _ur;

	std::vector<POINT*> _pnts;
	bool _is_leaf;
	/**
	 * maximum number of points per leaf
	 */
	static std::size_t _max_points_per_node;
};

template <typename POINT> std::size_t OctTree<POINT>::_max_points_per_node = 0;

} // end namespace GeoLib

#endif /* OCTTREE_H_ */
