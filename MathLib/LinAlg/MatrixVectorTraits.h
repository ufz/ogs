/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include "MatrixSpecifications.h"

namespace MathLib
{
template<typename Matrix>
struct MatrixVectorTraits;
}

#define SPECIALIZE_MATRIX_VECTOR_TRAITS(MATVEC, IDX) \
    template<> struct MatrixVectorTraits<MATVEC> { \
        using Index = IDX; \
        static std::unique_ptr<MATVEC> newInstance(); \
        static std::unique_ptr<MATVEC> newInstance(MATVEC const& A); \
        static std::unique_ptr<MATVEC> newInstance(MatrixSpecifications const& spec); \
    };


#ifdef USE_PETSC

#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"

namespace MathLib
{
SPECIALIZE_MATRIX_VECTOR_TRAITS(PETScMatrix, PETScMatrix::IndexType);
SPECIALIZE_MATRIX_VECTOR_TRAITS(PETScVector, PETScVector::IndexType);
}


#elif defined(OGS_USE_EIGEN)

#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"

namespace MathLib
{
SPECIALIZE_MATRIX_VECTOR_TRAITS(EigenMatrix, EigenMatrix::IndexType)
SPECIALIZE_MATRIX_VECTOR_TRAITS(EigenVector, EigenVector::IndexType)
}

#endif

#undef SPECIALIZE_MATRIX_VECTOR_TRAITS
