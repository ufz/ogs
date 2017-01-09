/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-22
 * \brief  Definition of the SimplePolygonTree class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
 * polygons. Every node of the SimplePolygonTree represents a polygon. A
 * child node c of a parent node p mean that the polygon represented by c
 * is contained in the polygon represented by p.
 */
class SimplePolygonTree
{
public:
    /** Creates a node of a tree containing a simple polygon.
     * @param polygon the polygon represented by this tree node
     * @param parent pointer to the parent node within the tree or nullptr
     * (if SimplePolygonTree node is the root node of the tree)
     */
    SimplePolygonTree(Polygon* polygon, SimplePolygonTree* parent);
    /** Destructor: Attention: does not destroy the polygon! */
    virtual ~SimplePolygonTree();
    /** Checks if the polygon represented by the given polygon tree node
     * is inside this node polygon.
     */
    bool isPolygonInside (const SimplePolygonTree* polygon_tree) const;
    /** Either insert the given SimplePolygonTree in one of the existing
     * children or as a new child.
     */
    void insertSimplePolygonTree (SimplePolygonTree* polygon_tree);

    /**
     * get the polygon represented by the tree node
     * @return the polygon
     */
    const Polygon* getPolygon () const;

    /** returns the number of children */
    std::size_t getNumberOfChildren() const { return _children.size(); }

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
    std::list<SimplePolygonTree*> _children;

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
    typedef typename std::list<POLYGONTREETYPE*>::const_iterator CIT;
    typedef typename std::list<POLYGONTREETYPE*>::iterator IT;
    for (CIT it0(list_of_simple_polygon_hierarchies.begin());
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
