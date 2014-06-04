/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-22
 * \brief  Definition of the SimplePolygonTree class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SIMPLEPOLYGONTREE_H_
#define SIMPLEPOLYGONTREE_H_

#include <list>

// GeoLib
#include "Polygon.h"

namespace GeoLib
{
/**
 * \brief This class computes and stores the topological relations between
 * polygons. Every node of the SimplePolygonTree represents a polygon.
 */
class SimplePolygonTree
{
public:
	SimplePolygonTree(Polygon* polygon, SimplePolygonTree* parent);
	virtual ~SimplePolygonTree();

	bool isPolygonInside (const SimplePolygonTree* polygon_tree) const;
	void insertSimplePolygonTree (SimplePolygonTree* polygon_tree);

	/**
	 * get the polygon represented by the tree node
	 * @return the polygon
	 */
	const Polygon* getPolygon () const;

	/** returns the number of childs */
	std::size_t getNChilds() const { return _childs.size(); }

protected:
	/**
	 * the polygon this node stands for
	 */
	Polygon* _node_polygon;
	/**
	 * the polygon represented by this node is contained in the
	 * polygon represented by the parent node in the tree
	 */
	SimplePolygonTree* _parent;
	/**
	 * list of polygons (represented by SimplePolygonTree nodes) contained
	 * in the _node_polygon
	 */
	std::list<SimplePolygonTree*> _childs;

private:
	void setParent(SimplePolygonTree* parent)
	{
		_parent = parent;
	}
};

/**
 * creates from a list of simple polygons a list of trees (SimplePolygonTrees)
 */
template <typename POLYGONTREETYPE>
void createPolygonTrees (std::list<POLYGONTREETYPE*>& list_of_simple_polygon_hierarchies)
{
	typedef typename std::list<POLYGONTREETYPE*>::const_iterator IT;
	for (IT it0(list_of_simple_polygon_hierarchies.begin());
		it0 != list_of_simple_polygon_hierarchies.end(); ++it0) {
		IT it1 = list_of_simple_polygon_hierarchies.begin();
		while (it1 != list_of_simple_polygon_hierarchies.end()) {
			if (it0 == it1) { // don't check same polygons
				++it1;
				// skip test if it1 points to the end after increment
				if (it1 == list_of_simple_polygon_hierarchies.end())
					break;
			}
			if ((*it0)->isPolygonInside(*it1)) {
				(*it0)->insertSimplePolygonTree(*it1);
				it1 = list_of_simple_polygon_hierarchies.erase(it1);
			} else {
				++it1;
			}
		}
	}
}


} // end namespace GeoLib

#endif /* SIMPLEPOLYGONTREE_H_ */
