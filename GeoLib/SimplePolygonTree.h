/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file SimplePolygonTree.h
 *
 *  Created on 2010-06-22 by Thomas Fischer
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
 *
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
 * @param list_of_simple_polygon_trees
 */
template <typename POLYGONTREETYPE>
void createPolygonTrees (std::list<POLYGONTREETYPE*>& list_of_simple_polygon_hierarchies)
{
	typename std::list<POLYGONTREETYPE*>::iterator it0 (list_of_simple_polygon_hierarchies.begin()), it1;
	while (it0 != list_of_simple_polygon_hierarchies.end()) {
		it1 = it0;
		it1++;
		while (it1 != list_of_simple_polygon_hierarchies.end()) {
			if ((*it0)->isPolygonInside(*it1)) {
				(*it0)->insertSimplePolygonTree(*it1);
				it1 = list_of_simple_polygon_hierarchies.erase(it1);
			} else {
				if ((*it1)->isPolygonInside(*it0)) {
					(*it1)->insertSimplePolygonTree(*it0);
					it0 = list_of_simple_polygon_hierarchies.erase(it0);
				}

				it1++;
			}
		}
		it0++;
	}
}


} // end namespace GeoLib

#endif /* SIMPLEPOLYGONTREE_H_ */
