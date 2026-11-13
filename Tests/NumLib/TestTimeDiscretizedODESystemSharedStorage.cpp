// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// Measures the heap memory taken by the global matrices and vector of
// TimeDiscretizedODESystem instances with and without shared storage.

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <optional>
#include <type_traits>
#include <vector>

#if defined(__has_feature)
#if __has_feature(address_sanitizer)
#define OGS_TEST_HEAP_REPLACED_BY_SANITIZER
#endif
#endif
#if defined(__SANITIZE_ADDRESS__)
#define OGS_TEST_HEAP_REPLACED_BY_SANITIZER
#endif

#if defined(__GLIBC__)
#include <malloc.h>
#endif

#include "BaseLib/Logging.h"
#include "NumLib/ODESolver/ODESystem.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"

// The sparsity pattern and the measured sizes are those of the Eigen matrices.
#ifndef USE_PETSC

namespace
{
//! An ODE that only provides matrix specifications with a sparsity pattern;
//! nothing is assembled.
//! Being a Newton ODE system, it is a Picard one, too.
class SparsePatternODE final
    : public NumLib::ODESystem<
          NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
          NumLib::NonlinearSolverTag::Newton>
{
public:
    SparsePatternODE(std::size_t const n_rows,
                     GlobalIndexType const n_nonzeros_per_row)
        : n_rows_(n_rows)
    {
        sparsity_pattern_.number_non_zeros_per_row.assign(n_rows,
                                                          n_nonzeros_per_row);
    }

    void preAssemble(double const /*t*/, double const /*dt*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& /*x*/,
                  std::vector<GlobalVector*> const& /*x_prev*/,
                  int const /*process_id*/, GlobalMatrix& /*M*/,
                  GlobalMatrix& /*K*/, GlobalVector& /*b*/) override
    {
    }

    void assembleWithJacobian(double const /*t*/, double const /*dt*/,
                              std::vector<GlobalVector*> const& /*x*/,
                              std::vector<GlobalVector*> const& /*x_prev*/,
                              int const /*process_id*/, GlobalVector& /*b*/,
                              GlobalMatrix& /*Jac*/) override
    {
    }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        int const /*process_id*/) const override
    {
        return {n_rows_, n_rows_, nullptr, &sparsity_pattern_};
    }

    bool isLinear() const override { return true; }

    bool requiresNormalization() const override { return false; }

private:
    std::size_t const n_rows_;
    GlobalSparsityPattern sparsity_pattern_;
};

//! Bytes of heap memory in use, including the chunks allocated by mmap, or
//! no value if that cannot be measured on this platform or the allocator is
//! replaced by the address sanitizer.
std::optional<std::size_t> heapBytesInUse()
{
#if defined(__GLIBC__) && !defined(OGS_TEST_HEAP_REPLACED_BY_SANITIZER) && \
    (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
    auto const info = mallinfo2();
    return info.uordblks + info.hblkhd;
#else
    return std::nullopt;
#endif
}

std::ptrdiff_t heapBytesSince(std::size_t const reference)
{
    return static_cast<std::ptrdiff_t>(*heapBytesInUse()) -
           static_cast<std::ptrdiff_t>(reference);
}

template <typename NLTagType>
class NumLibTimeDiscretizedODESystemSharedStorage : public ::testing::Test
{
public:
    static constexpr auto NLTag = NLTagType::value;
    using TDiscODESys = NumLib::TimeDiscretizedODESystem<
        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear, NLTag>;

    // Sizes of a scalar field with linear shape functions on a 3D hexahedral
    // mesh of about 100000 nodes.
    static constexpr std::size_t n_rows = 100'000;
    static constexpr GlobalIndexType n_nonzeros_per_row = 27;
    // Picard has M, K and b, Newton has the Jacobian and b.
    static constexpr std::size_t n_matrices =
        NLTag == NumLib::NonlinearSolverTag::Picard ? 2 : 1;

    void SetUp() override
    {
        if (!heapBytesInUse())
        {
            GTEST_SKIP() << "The heap usage cannot be measured on this "
                            "platform; glibc >= 2.33 without the address "
                            "sanitizer is required.";
        }
    }

    SparsePatternODE ode{n_rows, n_nonzeros_per_row};
    NumLib::BackwardEuler time_disc;
};

using NonlinearSolverTags = ::testing::Types<
    std::integral_constant<NumLib::NonlinearSolverTag,
                           NumLib::NonlinearSolverTag::Picard>,
    std::integral_constant<NumLib::NonlinearSolverTag,
                           NumLib::NonlinearSolverTag::Newton>>;

TYPED_TEST_SUITE(NumLibTimeDiscretizedODESystemSharedStorage,
                 NonlinearSolverTags);

TYPED_TEST(NumLibTimeDiscretizedODESystemSharedStorage, HeapMemoryUsage)
{
    using TDiscODESys = typename TestFixture::TDiscODESys;
    auto& ode = this->ode;
    auto& time_disc = this->time_disc;
    int const process_id = 0;
    auto const n_rows = TestFixture::n_rows;

    // Eigen's sparse matrix reserves a value and an index per nonzero, and
    // an outer index and a nonzero count per row.
    using StorageIndex = GlobalMatrix::RawMatrixType::StorageIndex;
    std::ptrdiff_t const n_nonzeros =
        n_rows * static_cast<std::size_t>(TestFixture::n_nonzeros_per_row);
    std::ptrdiff_t const matrix_bytes =
        n_nonzeros * (sizeof(double) + sizeof(StorageIndex)) +
        (2 * n_rows + 1) * sizeof(StorageIndex);
    std::ptrdiff_t const storage_bytes =
        TestFixture::n_matrices * matrix_bytes + n_rows * sizeof(double);
    // Covers the small allocations, e.g. of the shared pointers' control
    // blocks, and the freed chunks the allocator still counts as used.
    std::ptrdiff_t const tolerance = storage_bytes / 1000;

    // The first construction also makes some one-time allocations.
    {
        TDiscODESys const warm_up(process_id, ode, time_disc);
    }

    auto const reference = *heapBytesInUse();

    std::ptrdiff_t one_system = 0;
    {
        TDiscODESys const owner(process_id, ode, time_disc);
        one_system = heapBytesSince(reference);
    }
    EXPECT_LE(std::abs(one_system - storage_bytes), tolerance);
    EXPECT_LE(std::abs(heapBytesSince(reference)), tolerance);

    std::ptrdiff_t two_systems = 0;
    {
        TDiscODESys const first(process_id, ode, time_disc);
        TDiscODESys const second(process_id, ode, time_disc);
        two_systems = heapBytesSince(reference);
    }
    EXPECT_LE(std::abs(two_systems - 2 * storage_bytes), tolerance);
    EXPECT_LE(std::abs(heapBytesSince(reference)), tolerance);

    std::ptrdiff_t two_sharing_systems = 0;
    std::ptrdiff_t sharing_system_alone = 0;
    {
        auto owner = std::make_unique<TDiscODESys>(process_id, ode, time_disc);
        TDiscODESys const sharing(process_id, ode, time_disc, *owner);
        two_sharing_systems = heapBytesSince(reference);

        // The storage outlives the system that allocated it.
        owner.reset();
        sharing_system_alone = heapBytesSince(reference);
    }
    EXPECT_LE(std::abs(two_sharing_systems - storage_bytes), tolerance);
    EXPECT_LE(std::abs(sharing_system_alone - storage_bytes), tolerance);
    EXPECT_LE(std::abs(heapBytesSince(reference)), tolerance);

    INFO(
        "Heap memory of {:d} matrices and a vector of {:d} rows with {:d} "
        "nonzeros per row: expected {:d} B; one system {:d} B, two systems "
        "{:d} B, two sharing systems {:d} B.",
        TestFixture::n_matrices, n_rows, TestFixture::n_nonzeros_per_row,
        storage_bytes, one_system, two_systems, two_sharing_systems);
}
}  // namespace
#endif  // USE_PETSC
