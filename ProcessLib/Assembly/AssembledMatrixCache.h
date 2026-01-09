// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"

namespace ProcessLib
{
//! Stores assembled global M, K, b matrices/vectors for later reuse.
struct AssembledMatrixCache final
{
    //! Access the stored data.
    std::tuple<GlobalMatrix&, GlobalMatrix&, GlobalVector&> MKb() const
    {
        return std::tie(*M_, *K_, *b_);
    }

    //! Check if data has been stored.
    bool hasMKb() const { return M_ && K_ && b_; }

    //! Store data.
    void storeMKb(GlobalMatrix const& M,
                  GlobalMatrix const& K,
                  GlobalVector const& b)
    {
        M_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(M);
        K_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(K);
        b_ = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(b);
    }

private:
    std::unique_ptr<GlobalMatrix> M_;
    std::unique_ptr<GlobalMatrix> K_;
    std::unique_ptr<GlobalVector> b_;
};

/// Provides global M, K, b matrices/vectors either from an AssembledMatrixCache
/// or new instances if the cache does not hold a value.
struct MKbFromCacheOrFromThis
{
    /// Indicates if returned data came form a cache or from the fields of this
    /// MKbFromCacheOrFromThis instance.
    enum class From
    {
        Cache,
        This
    };

    /// Returns matrices and vectors and where they came from either from the
    /// passed AssembledMatrixCache or from the fields of this
    /// MKbFromCacheOrFromThis instance.
    std::tuple<From, GlobalMatrix&, GlobalMatrix&, GlobalVector&> MKb(
        AssembledMatrixCache const& cache,
        MathLib::MatrixSpecifications const& mat_spec)
    {
        if (cache.hasMKb())
        {
            return std::tuple_cat(std::tuple{From::Cache}, cache.MKb());
        }

        M_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(mat_spec);
        K_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(mat_spec);
        b_ = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(mat_spec);

        return std::tuple<From, GlobalMatrix&, GlobalMatrix&, GlobalVector&>{
            From::This, *M_, *K_, *b_};
    }

private:
    /// Matrices/vectors returned if data cannot be taken from a cache.
    /// These matrices/vectors will be allocated only if they are needed. I.e.,
    /// if the cache is populated in a timestep they might remain nullptr.
    std::unique_ptr<GlobalMatrix> M_;
    std::unique_ptr<GlobalMatrix> K_;
    std::unique_ptr<GlobalVector> b_;
};
}  // namespace ProcessLib
