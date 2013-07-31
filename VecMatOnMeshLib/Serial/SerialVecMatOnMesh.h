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

#ifndef VECMATONMESH_H_
#define VECMATONMESH_H_

#include <vector>

#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"

#include "../VecMeshItems/MeshComponentMap.h"
#include "ForEachMeshItem.h"

namespace VecMatOnMeshLib
{

/**
 * Non-parallel version using default LinAlg
 */
class SerialVectorMatrixBuilder
{
public:
    typedef MathLib::DenseVector<double> VectorType;
    typedef MathLib::GlobalDenseMatrix<double> MatrixType; //TODO this should be CRSMatrix
    template <class T_MESHITEM, class T_TASK>
    using ForEachType = ForEachMeshItem<T_MESHITEM, T_TASK>;

public:
    VectorType* createVector(const MeshComponentMap &dist_layout)
    {
        VectorType* vec = new VectorType(dist_layout.size());
        return vec;
    }

    MatrixType* createMatrix(const MeshComponentMap &dist_layout)
    {
        MatrixType* mat = new MatrixType(dist_layout.size(),dist_layout.size()); //TODO sparse structure
        return mat;
    }

};

} // VecMatOnMeshLib

#endif /* VECMATONMESH_H_ */
