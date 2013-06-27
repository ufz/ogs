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

#ifndef IFOREACHMESHITEM_H_
#define IFOREACHMESHITEM_H_

#include <vector>

namespace VecMatOnMeshLib
{


template <class T_MESHITEM, class T_TASK>
class IForEachMeshItem
{
public:
    virtual void operator()(const std::vector<T_MESHITEM*> &vec_items, T_TASK &task) = 0;

    virtual ~IForEachMeshItem() {};

};

}

#endif /* IFOREACHMESHITEM_H_ */
