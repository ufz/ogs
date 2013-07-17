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
        BY_COMPONENT,	///< Ordering data first by component type
        BY_LOCATION   	///< Ordering data first by spatial location
    };
};

}

#endif //NUMBERINGTYPE_H_
