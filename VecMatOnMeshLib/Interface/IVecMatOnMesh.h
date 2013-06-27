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

#ifndef IVECMATONMESH_H_
#define IVECMATONMESH_H_

namespace VecMatOnMeshLib
{

class VectorComposition;

/**
 * \brief Interface of VecMatOnMesh functions
 *
 *
 */
template <class T_VEC, class T_MAT>
class IVecMatOnMesh
{
public:
    virtual ~IVecMatOnMesh(){};

    virtual T_VEC* createVector(const VectorComposition &dist_layout) = 0;
    virtual T_MAT* createMatrix(const VectorComposition &dist_layout) = 0;
};

}

#endif /* IVECMATONMESH_H_ */
