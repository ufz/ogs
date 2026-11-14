// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// Measures the heap memory taken by the global matrices and vector of
// TimeDiscretizedODESystem instances with and without shared storage, and
// checks that sharing the storage does not change the solution of a staggered
// scheme.

#include <gtest/gtest.h>

#include <array>
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
#include "CreateTestLinearSolver.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/ODESolver/ConvergenceCriterionDeltaX.h"
#include "NumLib/ODESolver/FixedDampingStrategy.h"
#include "NumLib/ODESolver/NonlinearSolver.h"
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
/*! Two coupled nonlinear 1D diffusion-reaction equations, one per process id,
 * with the same sparsity pattern but different coefficients. The equation of
 * process p for x = x_p, coupled to y = x_{1-p}, reads
 *
 *   c_p dx/dt - D_p(y) Laplace(x) + k_p x^2 = s_p (1 + y),  D_p = 1 + a_p y^2,
 *
 * with a Dirichlet condition at the first node and a no-flux condition at the
 * last one.
 */
class CoupledNonlinearODE final
    : public NumLib::ODESystem<
          NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
          NumLib::NonlinearSolverTag::Newton>
{
public:
    static constexpr std::size_t n_rows = 20;

    CoupledNonlinearODE()
    {
        sparsity_pattern_.number_non_zeros_per_row.assign(n_rows, 3);
        for (int p = 0; p < 2; ++p)
        {
            known_solutions_[p].push_back({{0}, {1.0 + p}});
        }
    }

    void preAssemble(double const /*t*/, double const /*dt*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& x,
                  std::vector<GlobalVector*> const& /*x_prev*/,
                  int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        storage_addresses[process_id] = {&M, &b};
        auto const& [c, a, k, s] = coefficients_[process_id];
        auto const& u = *x[process_id];
        auto const& y = *x[1 - process_id];
        for (std::size_t i = 0; i < n_rows; ++i)
        {
            auto const row = static_cast<GlobalIndexType>(i);
            M.add(row, row, c);
            addDiffusion(K, i, 1. + a * y[i] * y[i]);
            K.add(row, row, k * u[i]);  // Picard linearisation of k x^2.
            b.add(row, s * (1. + y[i]));
        }
    }

    void assembleWithJacobian(double const /*t*/, double const dt,
                              std::vector<GlobalVector*> const& x,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        storage_addresses[process_id] = {&Jac, &b};
        auto const& [c, a, k, s] = coefficients_[process_id];
        auto const& u = *x[process_id];
        auto const& u_prev = *x_prev[process_id];
        auto const& y = *x[1 - process_id];
        for (std::size_t i = 0; i < n_rows; ++i)
        {
            auto const row = static_cast<GlobalIndexType>(i);
            double const D = 1. + a * y[i] * y[i];
            // b = -r, the negative residual.
            double r =
                c * (u[i] - u_prev[i]) / dt + k * u[i] * u[i] - s * (1. + y[i]);
            for (auto const j : neighbours(i))
            {
                r += D * (u[i] - u[j]);
            }
            b.add(row, -r);
            Jac.add(row, row, c / dt + 2. * k * u[i]);
            addDiffusion(Jac, i, D);
        }
    }

    std::vector<NumLib::IndexValueVector<Index>> const* getKnownSolutions(
        double const /*t*/, GlobalVector const& /*x*/,
        int const process_id) const override
    {
        return &known_solutions_[process_id];
    }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        int const /*process_id*/) const override
    {
        return {n_rows, n_rows, nullptr, &sparsity_pattern_};
    }

    bool isLinear() const override { return false; }

    bool requiresNormalization() const override { return false; }

    //! The matrix (M or Jacobian) and vector b the last assembly of each
    //! process id wrote to.
    std::array<std::pair<GlobalMatrix const*, GlobalVector const*>, 2>
        storage_addresses{};

private:
    static std::vector<std::size_t> neighbours(std::size_t const i)
    {
        std::vector<std::size_t> js;
        if (i > 0)
        {
            js.push_back(i - 1);
        }
        if (i + 1 < n_rows)
        {
            js.push_back(i + 1);
        }
        return js;
    }

    //! Adds D * (x_i - x_j) for the neighbours j of i to row i.
    static void addDiffusion(GlobalMatrix& A, std::size_t const i,
                             double const D)
    {
        auto const row = static_cast<GlobalIndexType>(i);
        for (auto const j : neighbours(i))
        {
            A.add(row, row, D);
            A.add(row, static_cast<GlobalIndexType>(j), -D);
        }
    }

    struct Coefficients
    {
        double c, a, k, s;
    };
    std::array<Coefficients, 2> const coefficients_{
        {{1.0, 0.5, 0.3, 1.0}, {2.0, 2.0, 1.5, 0.4}}};

    std::array<std::vector<NumLib::IndexValueVector<Index>>, 2>
        known_solutions_;
    GlobalSparsityPattern sparsity_pattern_;
};

template <NumLib::NonlinearSolverTag NLTag>
std::unique_ptr<NumLib::NonlinearSolver<NLTag>> createNonlinearSolver(
    GlobalLinearSolver& linear_solver)
{
    int const maxiter = 50;
    if constexpr (NLTag == NumLib::NonlinearSolverTag::Newton)
    {
        return std::make_unique<NumLib::NonlinearSolver<NLTag>>(
            linear_solver, maxiter,
            std::make_unique<NumLib::FixedDampingStrategy>(1.0));
    }
    else
    {
        return std::make_unique<NumLib::NonlinearSolver<NLTag>>(
            linear_solver, maxiter, /*anderson_depth=*/0, /*damping=*/1.0);
    }
}

/*! Solves the two equations of CoupledNonlinearODE with a staggered scheme,
 * as the time loop does: one nonlinear solver for both, the equations solved
 * one after the other in each coupling iteration.
 *
 * \return the solution of both equations at the end, and whether the two
 * equation systems assembled into the same matrix and vector.
 */
template <NumLib::NonlinearSolverTag NLTag>
std::pair<std::array<std::vector<double>, 2>, bool> solveStaggered(
    bool const share_storage)
{
    using TDiscODESys = NumLib::TimeDiscretizedODESystem<
        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear, NLTag>;
    std::size_t const n_rows = CoupledNonlinearODE::n_rows;

    CoupledNonlinearODE ode;
    std::array<NumLib::BackwardEuler, 2> time_discs;
    std::array<std::unique_ptr<TDiscODESys>, 2> systems;
    systems[0] = std::make_unique<TDiscODESys>(0, ode, time_discs[0]);
    systems[1] =
        share_storage
            ? std::make_unique<TDiscODESys>(1, ode, time_discs[1], *systems[0])
            : std::make_unique<TDiscODESys>(1, ode, time_discs[1]);

    auto linear_solver = createLinearSolver();
    auto nonlinear_solver = createNonlinearSolver<NLTag>(*linear_solver);
    std::array<std::unique_ptr<NumLib::ConvergenceCriterionDeltaX>, 2>
        convergence_criteria;
    for (auto& cc : convergence_criteria)
    {
        cc = std::make_unique<NumLib::ConvergenceCriterionDeltaX>(
            1e-12, std::nullopt, MathLib::VecNormType::NORM2);
    }

    std::array<GlobalVector, 2> x{GlobalVector(n_rows), GlobalVector(n_rows)};
    std::array<GlobalVector, 2> x_prev{GlobalVector(n_rows),
                                       GlobalVector(n_rows)};
    for (int p = 0; p < 2; ++p)
    {
        MathLib::LinAlg::set(x[p], 0.5 - 0.3 * p);
        MathLib::LinAlg::copy(x[p], x_prev[p]);
    }
    std::vector<GlobalVector*> const xs{&x[0], &x[1]};
    std::vector<GlobalVector*> const xs_prev{&x_prev[0], &x_prev[1]};

    double const dt = 0.1;
    for (auto& time_disc : time_discs)
    {
        time_disc.setInitialState(0.0);
    }
    for (int step = 1; step <= 3; ++step)
    {
        for (auto& time_disc : time_discs)
        {
            time_disc.nextTimestep(step * dt, dt);
        }
        for (int coupling_iteration = 0; coupling_iteration < 4;
             ++coupling_iteration)
        {
            for (int p = 0; p < 2; ++p)
            {
                nonlinear_solver->setEquationSystem(*systems[p],
                                                    *convergence_criteria[p]);
                auto const status =
                    nonlinear_solver->solve(xs, xs_prev, nullptr, p);
                EXPECT_TRUE(status.error_norms_met);
            }
        }
        for (int p = 0; p < 2; ++p)
        {
            MathLib::LinAlg::copy(x[p], x_prev[p]);
        }
    }

    std::array<std::vector<double>, 2> solution;
    for (int p = 0; p < 2; ++p)
    {
        for (std::size_t i = 0; i < n_rows; ++i)
        {
            solution[p].push_back(x[p][i]);
        }
    }
    return {solution, ode.storage_addresses[0] == ode.storage_addresses[1]};
}

template <typename NLTagType>
class NumLibTimeDiscretizedODESystemSharedStorageStaggered
    : public ::testing::Test
{
};

TYPED_TEST_SUITE(NumLibTimeDiscretizedODESystemSharedStorageStaggered,
                 NonlinearSolverTags);

// The solution must be bit-identical with and without shared storage, because
// each equation zeroes the shared matrices and vector before it assembles and
// has consumed them before the other one assembles.
TYPED_TEST(NumLibTimeDiscretizedODESystemSharedStorageStaggered,
           SameSolutionAsWithoutSharing)
{
    constexpr auto NLTag = TypeParam::value;

    auto const [own_storage_solution, own_storage_is_shared] =
        solveStaggered<NLTag>(false);
    auto const [shared_storage_solution, shared_storage_is_shared] =
        solveStaggered<NLTag>(true);

    EXPECT_FALSE(own_storage_is_shared);
    EXPECT_TRUE(shared_storage_is_shared);

    // Both equations have different, nontrivial solutions, otherwise stale
    // values of the other equation would go unnoticed.
    EXPECT_NE(own_storage_solution[0], own_storage_solution[1]);
    for (int p = 0; p < 2; ++p)
    {
        // Exact comparison on purpose.
        EXPECT_EQ(own_storage_solution[p], shared_storage_solution[p])
            << "process " << p;
    }
}
}  // namespace
#endif  // USE_PETSC
