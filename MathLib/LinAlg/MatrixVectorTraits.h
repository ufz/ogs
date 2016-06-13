/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_MATRIX_VECTOR_TRAITS_H
#define MATHLIB_MATRIX_VECTOR_TRAITS_H

#include<memory>

namespace MathLib
{
template<typename Matrix>
struct MatrixVectorTraits;

struct MatrixSpecifications;
}

#define SPECIALIZE_MATRIX_VECTOR_TRAITS(MATVEC, IDX) \
    template<> struct MatrixVectorTraits<MATVEC> { \
        using Index = IDX; \
        static std::unique_ptr<MATVEC> newInstance(); \
        static std::unique_ptr<MATVEC> newInstance(MATVEC const& A); \
        static std::unique_ptr<MATVEC> newInstance(MatrixSpecifications const& spec); \
    };


#ifdef OGS_USE_EIGEN

#include<Eigen/Core>

namespace MathLib
{
SPECIALIZE_MATRIX_VECTOR_TRAITS(Eigen::MatrixXd, Eigen::MatrixXd::Index)
SPECIALIZE_MATRIX_VECTOR_TRAITS(Eigen::VectorXd, Eigen::VectorXd::Index)
}

#endif


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

#endif // MATHLIB_MATRIX_VECTOR_TRAITS_H
