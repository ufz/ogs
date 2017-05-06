/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MatrixVectorTraits.h"
#include "MatrixSpecifications.h"

#ifdef USE_PETSC

namespace MathLib
{

std::unique_ptr<PETScMatrix>
MatrixVectorTraits<PETScMatrix>::
newInstance()
{
    return std::unique_ptr<PETScMatrix>(new PETScMatrix);
}

std::unique_ptr<PETScMatrix>
MatrixVectorTraits<PETScMatrix>::
newInstance(PETScMatrix const& A)
{
    return std::unique_ptr<PETScMatrix>(new PETScMatrix(A));
}

std::unique_ptr<PETScMatrix>
MatrixVectorTraits<PETScMatrix>::
newInstance(MatrixSpecifications const& spec)
{
    auto const nrows = spec.nrows;
    auto const ncols = spec.ncols;

    if (spec.sparsity_pattern)
    {
        // Assert that the misuse of the sparsity pattern is consistent.
        assert(spec.sparsity_pattern->size() == 1);

        auto const max_nonzeroes = spec.sparsity_pattern->front();

        PETScMatrixOption mat_opt;
        mat_opt.d_nz = max_nonzeroes;
        mat_opt.o_nz = max_nonzeroes;
        mat_opt.is_global_size = false;
        return std::unique_ptr<PETScMatrix>(
            new PETScMatrix(nrows, ncols, mat_opt));
    }
    else
        return std::unique_ptr<PETScMatrix>(
            new PETScMatrix(nrows, ncols));
}

std::unique_ptr<PETScVector>
MatrixVectorTraits<PETScVector>::
newInstance()
{
    return std::unique_ptr<PETScVector>(new PETScVector);
}

std::unique_ptr<PETScVector>
MatrixVectorTraits<PETScVector>::
newInstance(PETScVector const& x)
{
    return std::unique_ptr<PETScVector>(new PETScVector(x));
}

std::unique_ptr<PETScVector>
MatrixVectorTraits<PETScVector>::
newInstance(MatrixSpecifications const& spec)
{
    auto const is_global_size = false;

    if (spec.ghost_indices != nullptr) {
        return std::unique_ptr<PETScVector>(
            new PETScVector(spec.nrows, *spec.ghost_indices, is_global_size));
    } else {
        return std::unique_ptr<PETScVector>(
            new PETScVector(spec.nrows, is_global_size));
    }
}

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

namespace MathLib
{

std::unique_ptr<EigenMatrix>
MatrixVectorTraits<EigenMatrix>::
newInstance()
{
    return std::unique_ptr<EigenMatrix>(new EigenMatrix(0, 0)); // TODO default constructor
}

std::unique_ptr<EigenMatrix>
MatrixVectorTraits<EigenMatrix>::
newInstance(EigenMatrix const& A)
{
    return std::unique_ptr<EigenMatrix>(new EigenMatrix(A));
}

std::unique_ptr<EigenMatrix>
MatrixVectorTraits<EigenMatrix>::
newInstance(MatrixSpecifications const& spec)
{
    auto A = std::unique_ptr<EigenMatrix>(new EigenMatrix(spec.nrows));

    if (spec.sparsity_pattern)
        setMatrixSparsity(*A, *spec.sparsity_pattern);

    return A;
}

std::unique_ptr<EigenVector>
MatrixVectorTraits<EigenVector>::
newInstance()
{
    return std::unique_ptr<EigenVector>(new EigenVector);
}

std::unique_ptr<EigenVector>
MatrixVectorTraits<EigenVector>::
newInstance(EigenVector const& x)
{
    return std::unique_ptr<EigenVector>(new EigenVector(x));
}

std::unique_ptr<EigenVector>
MatrixVectorTraits<EigenVector>::
newInstance(MatrixSpecifications const& spec)
{
    return std::unique_ptr<EigenVector>(new EigenVector(spec.nrows));
}

} // namespace MathLib

#endif // defined(OGS_USE_EIGEN)
