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

#ifndef ITASK_H_
#define ITASK_H_

#include <cstddef>

namespace VecMatOnMeshLib
{

/**
 * Mesh item wise task interface
 *
 * \tparam T_MESH_ITEM    Mesh item type in MeshItemType.h
 */
template<class T_MESH_ITEM>
class ITask
{
public:
    virtual ~ITask() {};

    /**
     * do some task for the given mesh item
     *
     * @param item  Pointer to a mesh item
     * @param id    Item index
     */
    virtual void operator()(const T_MESH_ITEM* item, std::size_t id) = 0;
};

}

#endif /* ITASK_H_ */
