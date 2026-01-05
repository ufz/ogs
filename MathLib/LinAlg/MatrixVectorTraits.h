// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "MatrixSpecifications.h"

namespace MathLib
{
template <typename Matrix>
struct MatrixVectorTraits;
}

#define SPECIALIZE_MATRIX_VECTOR_TRAITS(MATVEC, IDX)                    \
    template <>                                                         \
    struct MatrixVectorTraits<MATVEC>                                   \
    {                                                                   \
        using Index = IDX;                                              \
        static std::unique_ptr<MATVEC> newInstance();                   \
        static std::unique_ptr<MATVEC> newInstance(MATVEC const& A);    \
        static std::unique_ptr<MATVEC> newInstance(                     \
            MatrixSpecifications const& spec);                          \
        static std::unique_ptr<MATVEC> newInstance(Index const length); \
    };

#ifdef USE_PETSC

#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"

namespace MathLib
{
SPECIALIZE_MATRIX_VECTOR_TRAITS(PETScMatrix, PETScMatrix::IndexType);
SPECIALIZE_MATRIX_VECTOR_TRAITS(PETScVector, PETScVector::IndexType);
}  // namespace MathLib

#else

#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"

namespace MathLib
{
SPECIALIZE_MATRIX_VECTOR_TRAITS(EigenMatrix, EigenMatrix::IndexType)
SPECIALIZE_MATRIX_VECTOR_TRAITS(EigenVector, EigenVector::IndexType)
}  // namespace MathLib

#endif

#undef SPECIALIZE_MATRIX_VECTOR_TRAITS
