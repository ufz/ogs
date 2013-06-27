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


#ifndef MESHITEM_H_
#define MESHITEM_H_

#include "MeshItemType.h"

namespace VecMatOnMeshLib
{

/**
 * Mesh item data
 *
 * This data are used in VectorComposition
 */
struct MeshItem
{
    std::size_t          mesh_id;
    MeshItemType::type   item_type;
    std::size_t          item_id;

    MeshItem(std::size_t meshid, MeshItemType::type itemtype, std::size_t itemid)
    : mesh_id(meshid), item_type(itemtype), item_id(itemid){};
};

inline bool operator<(const MeshItem& left, const MeshItem& right)
{
    if (left.mesh_id != right.mesh_id) return left.mesh_id < right.mesh_id;
    if (left.item_type != right.item_type) return left.item_type < right.item_type;
    return left.item_id < right.item_id;
}


}


#endif /* MESHITEM_H_ */
