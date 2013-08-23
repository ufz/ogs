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

#ifndef ASSEMBLERLIB_SERIALDENSEVECTORMATRIXBUILDER_H_
#define ASSEMBLERLIB_SERIALDENSEVECTORMATRIXBUILDER_H_

#include "VecMatOnMeshLib/Serial/SerialVectorMatrixBuilder.h"

#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"

namespace AssemblerLib
{

/// Serial vector/matrix builder using the GlobalDenseMatrix and the DenseVector
/// implementations.
///
/// \attention Using of GlobalDenseMatrix requires substantial amounts of
/// memory and this implementation is meant to be used for test purpose only.
typedef SerialVectorMatrixBuilder<
        MathLib::GlobalDenseMatrix<double>,
        MathLib::DenseVector<double>
    > SerialDenseVectorMatrixBuilder;

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALVECTOR_MATRIXBUILDER_H_
