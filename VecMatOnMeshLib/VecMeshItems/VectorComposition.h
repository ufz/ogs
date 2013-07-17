/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef DATAARRANGEMENT_H_
#define DATAARRANGEMENT_H_

#include <vector>

#include "MeshSubsets.h"
#include "OrderingType.h"
#include "MeshItem.h"
#include "MeshItemDataPositionDictionary.h"

namespace VecMatOnMeshLib
{

/**
 * VectorComposition class defines arrangement of distributed data over mesh items
 */
class VectorComposition
{
public:
	/**
	 *
	 * @param vec_comp_dis    a vector of component distributions
	 * The size of the vector means the number of components in the vector.
	 * @param ordering        type of ordering values in a vector
	 */
    VectorComposition(const std::vector<MeshSubsets*> &vec_comp_dis, OrderingType::type numbering);

    /// return the size of the vector
    std::size_t size() const { return _dict.size(); }

    /// return a vector of component IDs on a given mesh item
    std::vector<std::size_t> getComponentIDs(const Location &item) const;

    /// find a mesh item on which a given data position is assigned
    Location getMeshItem(unsigned dataID) const;

    /// find a data position from a given mesh item and component ID
    std::size_t getDataID(const Location &item, unsigned compID) const;

    /**
     * return a vector of data positions related to a given mesh item
     *
     * @param item   Location
     * @return a vector of data positions
     * If there is more than one component on a given mesh item, the function returns
     * a vector containing data positions for each component.
     */
    std::vector<std::size_t> getDataIDList(const Location &item) const;

    /**
     * return a vector of data positions corresponding to give items
     *
     * @param vec_items       a vector of mesh items
     * @param list_ordering   ordering type of data positions in a resulted vector
     * @return
     */
    std::vector<std::size_t> getDataIDList(const std::vector<Location> &vec_items, OrderingType::type list_ordering) const;

//for debugging
    const MeshitemDataPositionDictionary& getDictionary() const {return _dict; }
    void print();

private:
    void numberingByComponent(std::size_t offset=0);
    void numberingByMeshItems(std::size_t offset=0);

private:
    const OrderingType::type _ordering_type;
    MeshitemDataPositionDictionary _dict;
};

} // VecMatOnMeshLib

#endif /* DATAARRANGEMENT_H_ */
