/**
* \file OrderingType.h
* \author Norihiro Watanabe
* \date 2012-08-03
* \brief
*
* \copyright
* Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
* Distributed under a Modified BSD License.
* See accompanying file LICENSE.txt or
* http://www.opengeosys.org/project/license
*
*/

#ifndef NUMBERINGTYPE_H_
#define NUMBERINGTYPE_H_

namespace VecMatOnMeshLib
{

/**
* \brief Type of numbering distributed data index
*/
struct OrderingType
{
    enum type
    {
        BY_COMPONENT_TYPE, ///< Ordering data first by component type
        BY_MESH_ITEM_ID ///< Ordering data first by mesh item ID
    };
};

}

#endif //NUMBERINGTYPE_H_
