/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_SERIALVECTORMATRIXBUILDER_H_
#define ASSEMBLERLIB_SERIALVECTORMATRIXBUILDER_H_

#include "AssemblerLib/MeshComponentMap.h"

namespace AssemblerLib
{

template <typename MatrixType_, typename VectorType_>
class SerialVectorMatrixBuilder
{
public:
    typedef VectorType_ VectorType;
    typedef MatrixType_ MatrixType;

public:
    static
    VectorType* createVector(const MeshComponentMap &dist_layout)
    {
        VectorType* vec = new VectorType(dist_layout.size());
        return vec;
    }

    static
    MatrixType* createMatrix(const MeshComponentMap &dist_layout)
    {
        MatrixType* mat = new MatrixType(dist_layout.size());
        return mat;
    }

};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALVECTORMATRIXBUILDER_H_
