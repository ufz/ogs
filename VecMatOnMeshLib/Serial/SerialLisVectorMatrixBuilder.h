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

#ifndef ASSEMBLERLIB_SERIALLISVECTORMATRIXBUILDER_H_
#define ASSEMBLERLIB_SERIALLISVECTORMATRIXBUILDER_H_

#include "VecMatOnMeshLib/Serial/SerialVectorMatrixBuilder.h"

#include "MathLib/LinAlg/Lis/LisMatrix.h"
#include "MathLib/LinAlg/Lis/LisVector.h"

namespace AssemblerLib
{

/// Serial vector/matrix builder using LIS vectors and matrices.
typedef SerialVectorMatrixBuilder<
        MathLib::LisMatrix,
        MathLib::LisVector
    > SerialLisVectorMatrixBuilder;

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALLISVECTORMATRIXBUILDER_H_
