/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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
    /// Create matrix of given size. Any additional arguments are directly
    /// passed to the vector constructor.
    template <typename ...Args_>
    static
    VectorType* createVector(std::size_t const size, Args_&&... args)
    {
        return new VectorType(size, std::forward<Args_>(args)...);
    }

    /// Create a matrix of given size. Any additional arguments are directly
    /// passed to the matrix constructor.
    template <typename ...Args_>
    static
    MatrixType* createMatrix(std::size_t const size, Args_&&... args)
    {
        return new MatrixType(size, std::forward<Args_>(args)...);
    }

};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALVECTORMATRIXBUILDER_H_
