// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <range/v3/view/zip.hpp>

#include "BaseLib/Macros.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "ProcessLib/Assembly/MatrixAssemblyStats.h"

namespace ProcessLib::Assembly
{
namespace detail
{
#ifdef USE_PETSC
inline GlobalIndexType transformToNonGhostIndex(GlobalIndexType const i,
                                                GlobalIndexType const n)
{
    if (i == -n)
    {
        return 0;
    }
    else
    {
        return std::abs(i);
    }
}
#else
inline GlobalIndexType transformToNonGhostIndex(GlobalIndexType const i,
                                                GlobalIndexType const /*n*/)
{
    return i;
}
#endif
}  // namespace detail

template <std::size_t Dim>
struct MatrixElementCacheEntry
{
    MatrixElementCacheEntry() = default;
    MatrixElementCacheEntry(MatrixElementCacheEntry const&) = default;
    MatrixElementCacheEntry(MatrixElementCacheEntry&&) = default;

    MatrixElementCacheEntry(std::array<GlobalIndexType, Dim> const& indices,
                            double value)
        : indices{indices}, value{value}
    {
    }

    MatrixElementCacheEntry& operator=(MatrixElementCacheEntry const&) =
        default;
    MatrixElementCacheEntry& operator=(MatrixElementCacheEntry&&) = default;

    std::array<GlobalIndexType, Dim> indices;
    double value;
};

template <std::size_t Dim>
class MatrixElementCache final
{
    static_assert(Dim == 1 || Dim == 2);
    using MatOrVec = std::conditional_t<Dim == 1, GlobalVector, GlobalMatrix>;

    static constexpr std::size_t cache_capacity = 1'000'000;

public:
    MatrixElementCache(MatOrVec& mat_or_vec, Stats& stats,
                       const int num_threads)
        : mat_or_vec_(mat_or_vec), stats_(stats), num_threads_(num_threads)
    {
        cache_.reserve(cache_capacity);
    }

    void add(std::vector<double> const& local_data,
             std::vector<GlobalIndexType> const& indices)
    {
        if (local_data.empty())
        {
            return;
        }

        if (num_threads_ == 1)
        {
            if constexpr (Dim == 2)
            {
                auto const num_r_c = indices.size();
                mat_or_vec_.add(
                    MathLib::RowColumnIndices<GlobalIndexType>{indices,
                                                               indices},
                    MathLib::toMatrix(local_data, num_r_c, num_r_c));
            }
            else
            {
                mat_or_vec_.add(indices, local_data);
            }
            return;
        }

        addToCache(local_data, indices);
    }

    ~MatrixElementCache() { addCacheToGlobal(); }

private:
    void addToCache(std::vector<double> const& values,
                    std::vector<GlobalIndexType> const& indices)
    {
        ensureEnoughSpace(values.size());

        addToCacheImpl(values, indices,
                       std::integral_constant<std::size_t, Dim>{});
    }

    // Overload for vectors.
    void addToCacheImpl(std::vector<double> const& values,
                        std::vector<GlobalIndexType> const& indices,
                        std::integral_constant<std::size_t, 1>)
    {
        auto const num_r_c = indices.size();

        for (std::size_t r_local = 0; r_local < num_r_c; ++r_local)
        {
            ++stats_.count;
            auto const value = values[r_local];

            if (value == 0)
            {
                continue;
            }
            else
            {
                ++stats_.count_nonzero;
            }

            auto const r_global = indices[r_local];
            cache_.emplace_back(std::array{r_global}, value);
        }
    }

    // Overload for matrices.
    void addToCacheImpl(std::vector<double> const& values,
                        std::vector<GlobalIndexType> const& indices,
                        std::integral_constant<std::size_t, 2>)
    {
        auto const num_r_c = indices.size();

        // Note: There is an implicit storage order assumption, here!
        auto const local_mat = MathLib::toMatrix(values, num_r_c, num_r_c);

        for (std::size_t r_local = 0; r_local < num_r_c; ++r_local)
        {
            auto const r_global = indices[r_local];

            for (std::size_t c_local = 0; c_local < num_r_c; ++c_local)
            {
                ++stats_.count;
                auto const value = local_mat(r_local, c_local);

                // TODO skipping zero values sometimes does not work
                // together with the Eigen SparseLU linear solver. See also
                // https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/4556#note_125561
#if 0
                if (value == 0)
                {
                    continue;
                }
#endif
                ++stats_.count_nonzero;

                auto const c_global = indices[c_local];
                cache_.emplace_back(std::array{r_global, c_global}, value);
            }
        }
    }

    void ensureEnoughSpace(std::size_t const space_needed)
    {
        // Fast path - inline this check
        if (cache_.size() + space_needed <= cache_.capacity()) [[likely]]
        {
            return;
        }

        // Slow path - keep in separate function to avoid code bloat
        ensureEnoughSpaceSlow(space_needed);
    }

    OGS_NO_INLINE void ensureEnoughSpaceSlow(std::size_t const space_needed)
    {
        if (cache_.capacity() < cache_capacity)
        {
            cache_.reserve(cache_capacity);
        }

        auto const size_initial = cache_.size();
        auto const cap_initial = cache_.capacity();

        if (size_initial + space_needed <= cap_initial)
        {
            return;
        }

        addCacheToGlobal();

        // ensure again that there is enough capacity (corner case, if initial
        // capacity is too small because of whatever arcane reason)
        auto const size_after = cache_.size();
        auto const cap_after = cache_.capacity();

        if (size_after + space_needed > cap_after)
        {
            cache_.reserve(size_after + 2 * space_needed);
        }
    }

    void addCacheToGlobal()
    {
        if (cache_.empty())
        {
            return;
        }

#pragma omp critical
        {
            if constexpr (Dim == 2)
            {
                auto const n_cols = mat_or_vec_.getNumberOfColumns();

                // TODO would be more efficient if our global matrix and vector
                // implementations supported batch addition of matrix elements
                // with arbitrary indices (not restricted to (n x m) shaped
                // submatrices).
                for (auto const [rc, value] : cache_)
                {
                    auto const [r, c] = rc;

                    auto const c_no_ghost =
                        detail::transformToNonGhostIndex(c, n_cols);

                    mat_or_vec_.add(r, c_no_ghost, value);
                }
            }
            else
            {
                // TODO batch addition would be more efficient. That needs the
                // refactoring of the matrix element cache.
                for (auto const [r, value] : cache_)
                {
                    mat_or_vec_.add(r.front(), value);
                }
            }
        }

        stats_.count_global += cache_.size();
        cache_.clear();
    }

    std::vector<MatrixElementCacheEntry<Dim>> cache_;
    MatOrVec& mat_or_vec_;
    Stats& stats_;
    int num_threads_;
};

template <int NumMatrices>
class MultiMatrixElementCache;

// Specialization with 1 matrix ("Jac")
template <>
class MultiMatrixElementCache<1> final
{
public:
    MultiMatrixElementCache(GlobalVector& b, GlobalMatrix& Jac,
                            MultiStats<1>& stats, const int num_threads)
        : cache_b_(b, stats.b, num_threads),
          cache_Jac_(Jac, stats.Jac, num_threads)
    {
    }

    void add(std::vector<double> const& local_b_data,
             std::vector<double> const& local_Jac_data,
             std::vector<GlobalIndexType> const& indices)
    {
        cache_b_.add(local_b_data, indices);
        cache_Jac_.add(local_Jac_data, indices);
    }

private:
    MatrixElementCache<1> cache_b_;
    MatrixElementCache<2> cache_Jac_;
};

// Specialization with 2 matrices ("M", "K")
template <>
class MultiMatrixElementCache<2> final
{
public:
    MultiMatrixElementCache(GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
                            MultiStats<2>& stats, const int num_threads)
        : cache_M_(M, stats.M, num_threads),
          cache_K_(K, stats.K, num_threads),
          cache_b_(b, stats.b, num_threads)
    {
    }

    void add(std::vector<double> const& local_M_data,
             std::vector<double> const& local_K_data,
             std::vector<double> const& local_b_data,
             std::vector<GlobalIndexType> const& indices)
    {
        cache_M_.add(local_M_data, indices);
        cache_K_.add(local_K_data, indices);
        cache_b_.add(local_b_data, indices);
    }

private:
    MatrixElementCache<2> cache_M_;
    MatrixElementCache<2> cache_K_;
    MatrixElementCache<1> cache_b_;
};
}  // namespace ProcessLib::Assembly
