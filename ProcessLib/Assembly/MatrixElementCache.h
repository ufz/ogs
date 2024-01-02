/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>

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
class ConcurrentMatrixView
{
    static_assert(Dim == 1 || Dim == 2);
    using GlobalMatOrVec =
        std::conditional_t<Dim == 1, GlobalVector, GlobalMatrix>;

public:
    explicit ConcurrentMatrixView(GlobalMatOrVec& mat_or_vec)
        : mat_or_vec_{mat_or_vec}
    {
    }

    void add(std::vector<MatrixElementCacheEntry<Dim>> const& entries)
    {
        std::lock_guard<std::mutex> const lock(mutex_);

        if constexpr (Dim == 2)
        {
            auto const n_cols = mat_or_vec_.getNumberOfColumns();

            // TODO would be more efficient if our global matrix and vector
            // implementations supported batch addition of matrix elements with
            // arbitrary indices (not restricted to (n x m) shaped submatrices).
            for (auto const [rc, value] : entries)
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
            for (auto const [r, value] : entries)
            {
                mat_or_vec_.add(r.front(), value);
            }
        }
    }

private:
    std::mutex mutex_;
    GlobalMatOrVec& mat_or_vec_;
};

explicit ConcurrentMatrixView(GlobalVector&) -> ConcurrentMatrixView<1>;
explicit ConcurrentMatrixView(GlobalMatrix&) -> ConcurrentMatrixView<2>;

template <std::size_t Dim>
class MatrixElementCache final
{
    static_assert(Dim == 1 || Dim == 2);
    using GlobalMatView = ConcurrentMatrixView<Dim>;

    static constexpr std::size_t cache_capacity = 1'000'000;

public:
    MatrixElementCache(GlobalMatView& mat_or_vec, Stats& stats)
        : mat_or_vec_(mat_or_vec), stats_(stats)
    {
        cache_.reserve(cache_capacity);
    }

    void add(std::vector<double> const& local_data,
             std::vector<GlobalIndexType> const& indices)
    {
        addToCache(local_data, indices);
    }

    ~MatrixElementCache() { addToGlobal(); }

private:
    void addToCache(std::vector<double> const& values,
                    std::vector<GlobalIndexType> const& indices)
    {
        if (values.empty())
        {
            return;
        }

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

                // TODO skipping zero values sometimes does not work together
                // with the Eigen SparseLU linear solver. See also
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
        auto const size_initial = cache_.size();
        auto const cap_initial = cache_.capacity();

        if (size_initial + space_needed <= cap_initial)
        {
            return;
        }

        addToGlobal();

        // ensure again that there is enough capacity (corner case, if initial
        // capacity is too small because of whatever arcane reason)
        auto const size_after = cache_.size();
        auto const cap_after = cache_.capacity();

        if (size_after + space_needed > cap_after)
        {
            cache_.reserve(size_after + 2 * space_needed);
        }
    }

    void addToGlobal()
    {
        mat_or_vec_.add(cache_);
        stats_.count_global += cache_.size();
        cache_.clear();
    }

    std::vector<MatrixElementCacheEntry<Dim>> cache_;
    GlobalMatView& mat_or_vec_;
    Stats& stats_;
};

class MultiMatrixElementCache final
{
    using GlobalMatrixView = ConcurrentMatrixView<2>;
    using GlobalVectorView = ConcurrentMatrixView<1>;

public:
    MultiMatrixElementCache(GlobalMatrixView& M, GlobalMatrixView& K,
                            GlobalVectorView& b, GlobalMatrixView& Jac,
                            MultiStats& stats)
        : cache_M_(M, stats.M),
          cache_K_(K, stats.K),
          cache_b_(b, stats.b),
          cache_Jac_(Jac, stats.Jac)
    {
    }

    void add(std::vector<double> const& local_M_data,
             std::vector<double> const& local_K_data,
             std::vector<double> const& local_b_data,
             std::vector<double> const& local_Jac_data,
             std::vector<GlobalIndexType> const& indices)
    {
        cache_M_.add(local_M_data, indices);
        cache_K_.add(local_K_data, indices);
        cache_b_.add(local_b_data, indices);
        cache_Jac_.add(local_Jac_data, indices);
    }

private:
    MatrixElementCache<2> cache_M_;
    MatrixElementCache<2> cache_K_;
    MatrixElementCache<1> cache_b_;
    MatrixElementCache<2> cache_Jac_;
};
}  // namespace ProcessLib::Assembly
