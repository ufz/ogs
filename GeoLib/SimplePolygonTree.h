/**
 * \file SimplePolygonTree.h
 *
 *  Created on 2010-06-22 by Thomas Fischer
 */

#ifndef SIMPLEPOLYGONTREE_H_
#define SIMPLEPOLYGONTREE_H_

#include "Polygon.h"
// FileIO

namespace GeoLib {

/**
 * \brief This class computes and stores the topological relations between
 * polygons and geometric objects like Point and Polyline.
 *
 * It is used to generate a proper input file for gmsh.
 */
class SimplePolygonTree {
public:
	SimplePolygonTree(const Polygon* polygon, SimplePolygonTree* parent = NULL);
	virtual ~SimplePolygonTree();

	const Polygon* getPolygon () const;
	const std::list<SimplePolygonTree*>& getChilds() const;
	const std::list<GeoObject*>& getGeoObjects () const;
	size_t getLevel () const;

	bool isPolygonInside (const SimplePolygonTree* polygon_tree) const;
	void insertSimplePolygonTree (SimplePolygonTree* polygon_tree);
//	void visitAndProcessNodes (FileIO::GMSHInterface& gmsh_io);

	bool isGeoObjInside (const GeoObject* geo_obj) const;
	void insertGeoObj (const GeoObject* geo_obj);

private:
	bool isPolylineInside (const Polyline* ply) const;
//	void _visitAndProcessNodes (FileIO::GMSHInterface& gmsh_io);
	const Polygon* _node;
	SimplePolygonTree* _parent;
	std::list<SimplePolygonTree*> _childs;
	std::list<const GeoObject*> _geo_objs;
};

/**
 * creates from a list of simple polygons a list
 * @param list_of_simple_polygon_trees
 */
void createPolygonTree (std::list<SimplePolygonTree*>& list_of_simple_polygon_trees);

}

#endif /* SIMPLEPOLYGONTREE_H_ */
