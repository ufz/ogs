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

#include <vector>

#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#include "MathLib/LinAlg/Lis/LisMatrix.h"

#include "AssemblerLib/MeshComponentMap.h"

namespace AssemblerLib
{

using AssemblerLib::MeshComponentMap;

/**
 * Non-parallel version using default LinAlg
 */
template <typename MatrixType_, typename VectorType_>
class SerialVectorMatrixBuilder
{
public:
    typedef VectorType_ VectorType;
    typedef MatrixType_ MatrixType;

    /// Executes task(item, index) for each item in input vector.
    /// Return values of the task call are ignored.
    ///
    /// \param vec_items   a vector of mesh item pointers
    /// \param task        a function that accepts a mesh item pointer
    ///                    and an index as arguments
    template <class T_MESHITEM, class T_TASK>
    static
    void
    forEachMeshItem(const std::vector<T_MESHITEM*> &vec_items, T_TASK &task)
    {
        for (std::size_t i=0; i<vec_items.size(); i++)
            task(vec_items[i], i);
    }

public:
    VectorType* createVector(const MeshComponentMap &dist_layout)
    {
        VectorType* vec = new VectorType(dist_layout.size());
        return vec;
    }

    MatrixType* createMatrix(const MeshComponentMap &dist_layout)
    {
        MatrixType* mat = new MatrixType(dist_layout.size());
        return mat;
    }

};

typedef SerialVectorMatrixBuilder<
        MathLib::GlobalDenseMatrix<double>,
        MathLib::DenseVector<double>
    > SerialDenseVectorMatrixBuilder;

typedef SerialVectorMatrixBuilder<
        MathLib::LisMatrix,
        MathLib::LisVector
    > SerialLisVectorMatrixBuilder;

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALVECTORMATRIXBUILDER_H_
