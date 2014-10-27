/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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
    VectorType* createVector(std::size_t const size)
    {
        VectorType* vec = new VectorType(size);
        return vec;
    }

    static
    MatrixType* createMatrix(std::size_t const size)
    {
        MatrixType* mat = new MatrixType(size);
        return mat;
    }

};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALVECTORMATRIXBUILDER_H_
