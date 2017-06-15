/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <random>
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"

template <typename Vector>
void fillVectorRandomly(Vector& x)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    std::uniform_real_distribution<double> rnd;

    using Index = typename MathLib::MatrixVectorTraits<Vector>::Index;
    Index const size = x.size();

    for (Index i = 0; i < size; ++i) {
        MathLib::setVector(x, i, rnd(random_number_generator));
    }
#ifdef USE_PETSC
    finalizeVectorAssembly(x);
#endif
}
