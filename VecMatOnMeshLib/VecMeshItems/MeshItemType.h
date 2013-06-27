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


#ifndef MESHITEMTYPE_H_
#define MESHITEMTYPE_H_

namespace VecMatOnMeshLib
{

struct MeshItemType
{
    enum type
    {
        Node,
        Edge,
        Face,
        Cell
    };
};

}

#endif /* MESHITEMTYPE_H_ */
