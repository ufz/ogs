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

#ifndef FOREACHMESHITEM_H_
#define FOREACHMESHITEM_H_

#include <vector>
#include <algorithm>

#include "../Interface/IForEachMeshItem.h"

namespace VecMatOnMeshLib
{

/**
 * Apply a function to given items
 */
template <class T_MESHITEM, class T_TASK>
class ForEachMeshItem : public IForEachMeshItem<T_MESHITEM, T_TASK>
{
public:
    /**
     * do some work for each item in a give list
     *
     * @param vec_items   a vector of mesh items
     * @param task        a function that accepts an element in the range and index as arguments
     * Its return value, if any, is ignored.
     */
    virtual void operator()(const std::vector<T_MESHITEM*> &vec_items, T_TASK &task)
    {
        for (std::size_t i=0; i<vec_items.size(); i++)
            task(vec_items[i], i);
    }

    virtual ~ForEachMeshItem() {};

};

} // VecMatOnMeshLib

#endif /* FOREACHMESHITEM_H_ */
