/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FluxCorrectedTransport.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <ranges>

#ifdef USE_PETSC
#include <petscerror.h>
#include <petscmat.h>
#include <petscvec.h>

#include "BaseLib/MPI.h"
#endif

namespace NumLib
{
namespace detail
{

template <typename MatrixVectorType>
std::unique_ptr<MatrixVectorType> newZeroedInstance(
    MathLib::MatrixSpecifications const& matrix_specification)
{
    auto result = MathLib::MatrixVectorTraits<MatrixVectorType>::newInstance(
        matrix_specification);
    result->setZero();
    return result;
}

#ifndef USE_PETSC
void calculateFluxCorrectedTransportSerial(
    [[maybe_unused]] const double t, const double dt,
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    const MathLib::MatrixSpecifications& matrix_specification, GlobalMatrix& M,
    GlobalMatrix& K, GlobalVector& b)
{
    auto D = newZeroedInstance<GlobalMatrix>(matrix_specification);
    auto F = newZeroedInstance<GlobalMatrix>(matrix_specification);

    // compute artificial diffusion operator D
    using RawMatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    auto& K_raw = K.getRawMatrix();
    K_raw *= -1;
    for (int k = 0; k < K_raw.outerSize(); ++k)
    {
        for (RawMatrixType::InnerIterator it(K_raw, k); it; ++it)
        {
            if (it.row() != it.col())
            {
                double const kij = it.value();
                double const kji = K.get(it.col(), it.row());
                double const dij = std::max({-kij, 0., -kji});
                D->setValue(it.row(), it.col(), dij);
            }
        }
    }
    auto& D_raw = D->getRawMatrix();
    D_raw -= (D_raw * Eigen::VectorXd::Ones(D_raw.cols())).eval().asDiagonal();

    // compute F
    for (int k = 0; k < D_raw.outerSize(); ++k)
    {
        for (RawMatrixType::InnerIterator it(D_raw, k); it; ++it)
        {
            double const dij = it.value();
            double const xj = x[process_id]->get(it.col());
            double const xi = x[process_id]->get(it.row());
            double const fij = -dij * (xj - xi);
            F->setValue(it.row(), it.col(), fij);
        }
    }

    auto& M_raw = M.getRawMatrix();
    for (int k = 0; k < M_raw.outerSize(); ++k)
    {
        for (RawMatrixType::InnerIterator it(M_raw, k); it; ++it)
        {
            double const mij = it.value();
            double const xdotj = (x[process_id]->get(it.col()) -
                                  x_prev[process_id]->get(it.col())) /
                                 dt;
            double const xdoti = (x[process_id]->get(it.row()) -
                                  x_prev[process_id]->get(it.row())) /
                                 dt;
            double const fij = -mij * (xdotj - xdoti);
            F->add(it.row(), it.col(), fij);
        }
    }

    auto P_plus = newZeroedInstance<GlobalVector>(matrix_specification);
    auto P_minus = newZeroedInstance<GlobalVector>(matrix_specification);
    auto Q_plus = newZeroedInstance<GlobalVector>(matrix_specification);
    auto Q_minus = newZeroedInstance<GlobalVector>(matrix_specification);
    auto R_plus = newZeroedInstance<GlobalVector>(matrix_specification);
    auto R_minus = newZeroedInstance<GlobalVector>(matrix_specification);

    auto& F_raw = F->getRawMatrix();
    for (int k = 0; k < F_raw.outerSize(); ++k)
    {
        for (RawMatrixType::InnerIterator it(F_raw, k); it; ++it)
        {
            if (it.row() != it.col())
            {
                double const fij = it.value();
                P_plus->add(it.row(), std::max(0., fij));
                P_minus->add(it.row(), std::min(0., fij));

                double const x_prev_i = x_prev[process_id]->get(it.row());
                double const x_prev_j = x_prev[process_id]->get(it.col());

                double const Q_plus_i = Q_plus->get(it.row());
                double Q_plus_i_tmp =
                    std::max({Q_plus_i, 0., x_prev_j - x_prev_i});
                Q_plus->set(it.row(), Q_plus_i_tmp);

                double const Q_minus_i = Q_minus->get(it.row());
                double Q_minus_i_tmp =
                    std::min({Q_minus_i, 0., x_prev_j - x_prev_i});
                Q_minus->set(it.row(), Q_minus_i_tmp);
            }
        }
    }

    Eigen::VectorXd const M_L =
        (M_raw * Eigen::VectorXd::Ones(M_raw.cols())).eval();
    for (auto k = R_plus->getRangeBegin(); k < R_plus->getRangeEnd(); ++k)
    {
        double const P_plus_i = P_plus->get(k);

        if (P_plus_i == 0.0)
        {
            continue;
        }

        double const mi = M_L(k);
        double const Q_plus_i = Q_plus->get(k);
        R_plus->set(k, std::min(1.0, mi * Q_plus_i / dt / P_plus_i));
    }

    for (auto k = R_minus->getRangeBegin(); k < R_minus->getRangeEnd(); ++k)
    {
        double const P_minus_i = P_minus->get(k);
        if (P_minus_i == 0.0)
        {
            continue;
        }
        double const mi = M_L(k);
        double const Q_minus_i = Q_minus->get(k);
        R_minus->set(k, std::min(1.0, mi * Q_minus_i / dt / P_minus_i));
    }

    for (int k = 0; k < F_raw.outerSize(); ++k)
    {
        for (RawMatrixType::InnerIterator it(F_raw, k); it; ++it)
        {
            double const fij = it.value();
            if (fij > 0.)
            {
                double const R_plus_i = R_plus->get(it.row());
                double const R_minus_j = R_minus->get(it.col());
                double const alpha_ij = std::min(R_plus_i, R_minus_j);

                b.add(it.row(), alpha_ij * fij);
            }
            else
            {
                double const R_minus_i = R_minus->get(it.row());
                double const R_plus_j = R_plus->get(it.col());
                double const alpha_ij = std::min(R_minus_i, R_plus_j);

                b.add(it.row(), alpha_ij * fij);
            }
        }
    }

    // compute low-order operator
    K_raw += D_raw;
    K_raw *= -1;

    // overwrite with the lumped mass matrix
    M.setZero();
    for (int k = 0; k < M.getNumberOfRows(); ++k)
    {
        M.setValue(k, k, M_L(k));
    }
}
#endif  // end of ifndef USE_PETSC

#ifdef USE_PETSC
void finalize(Mat& M)
{
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatAssemblyBegin(M, MatAssemblyType::MAT_FINAL_ASSEMBLY));
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatAssemblyEnd(M, MatAssemblyType::MAT_FINAL_ASSEMBLY));
}

void calculateFluxCorrectedTransportPETSc(
    [[maybe_unused]] const double t, const double dt,
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    const MathLib::MatrixSpecifications& matrix_specification, GlobalMatrix& M,
    GlobalMatrix& K, GlobalVector& b)
{
    BaseLib::MPI::Mpi mpi(PETSC_COMM_WORLD);
    if (mpi.size > 1)
    {
        OGS_FATAL(
            "The current parallel implementation of flux corrected transport "
            "supports only one MPI rank, but the simulation uses {} ranks",
            mpi.size);
    }
    Mat K_raw = K.getRawMatrix();
    Mat D;
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatDuplicate(K_raw, MAT_DO_NOT_COPY_VALUES, &D));
    PetscCallAbort(PETSC_COMM_WORLD, MatZeroEntries(D));
    finalize(D);
    PetscInt start, end;
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatGetOwnershipRange(K.getRawMatrix(), &start, &end));

    PetscCallAbort(PETSC_COMM_WORLD, MatScale(K_raw, -1.0));
    K.finalizeAssembly();
    // compute artificial diffusion operator D
    for (PetscInt i = start; i < end; ++i)
    {
        PetscInt number_of_columns;
        PetscInt const* columns;
        PetscScalar const* values;

        PetscCallAbort(PETSC_COMM_WORLD, MatGetRow(K_raw, i, &number_of_columns,
                                                   &columns, &values));

        auto const columns_span =
            std::span<PetscInt const>(columns, number_of_columns);
        auto const values_span =
            std::span<PetscScalar const>(values, number_of_columns);
        for (auto const [j, kij] : std::views::zip(columns_span, values_span))
        {
            if (i == j)  // skip diagonal entries
            {
                continue;
            }
            PetscScalar kji;
            // MatGetValue is not supported for all PETSc matrix types
            PetscCallAbort(PETSC_COMM_WORLD, MatGetValue(K_raw, j, i, &kji));
            PetscScalar dij = std::max({-kij, 0.0, -kji});
            PetscCallAbort(PETSC_COMM_WORLD,
                           MatSetValue(D, i, j, dij, INSERT_VALUES));
        }
        // MatRestoreRow() is required after MatGetRow()
        PetscCallAbort(
            PETSC_COMM_WORLD,
            MatRestoreRow(K_raw, i, &number_of_columns, &columns, &values));
    }
    finalize(D);
    auto row_sums = newZeroedInstance<GlobalVector>(matrix_specification);
    Vec row_sums_raw = row_sums->getRawVector();
    PetscCallAbort(PETSC_COMM_WORLD, MatGetRowSum(D, row_sums_raw));
    PetscCallAbort(PETSC_COMM_WORLD, VecScale(row_sums_raw, -1.0));
    PetscCallAbort(PETSC_COMM_WORLD, VecAssemblyBegin(row_sums_raw));
    PetscCallAbort(PETSC_COMM_WORLD, VecAssemblyEnd(row_sums_raw));
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatDiagonalSet(D, row_sums_raw, ADD_VALUES));
    finalize(D);
    // need to access x and x_prev in the following loops
    x[process_id]->setLocalAccessibleVector();
    x_prev[process_id]->setLocalAccessibleVector();

    // compute F
    Mat F;
    PetscCallAbort(PETSC_COMM_WORLD,
                   MatDuplicate(K.getRawMatrix(), MAT_DO_NOT_COPY_VALUES, &F));
    PetscCallAbort(PETSC_COMM_WORLD, MatZeroEntries(F));
    for (PetscInt i = start; i < end; ++i)
    {
        PetscInt cols;
        const PetscInt* col_indices;
        const PetscScalar* values;

        PetscCallAbort(PETSC_COMM_WORLD,
                       MatGetRow(D, i, &cols, &col_indices, &values));

        auto const columns = std::span<PetscInt const>(col_indices, cols);
        auto const values_D = std::span<PetscScalar const>(values, cols);

        for (auto const [j, dij] : std::views::zip(columns, values_D))
        {
            PetscScalar const xi = x[process_id]->get(i);
            PetscScalar const xj = x[process_id]->get(j);
            PetscScalar const fij = -dij * (xj - xi);
            PetscCallAbort(PETSC_COMM_WORLD,
                           MatSetValue(F, i, j, fij, INSERT_VALUES));
        }
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatRestoreRow(D, i, &cols, &col_indices, &values));
    }
    finalize(F);
    auto& M_raw = M.getRawMatrix();

    auto xdot = newZeroedInstance<GlobalVector>(matrix_specification);
    xdot->setLocalAccessibleVector();
    for (PetscInt i = start; i < end; ++i)
    {
        xdot->add(i, (x[process_id]->get(i) - x_prev[process_id]->get(i)) / dt);
    }
    xdot->finalizeAssembly();
    xdot->setLocalAccessibleVector();

    for (PetscInt i = start; i < end; ++i)
    {
        PetscInt cols;
        PetscInt const* col_indices;
        PetscScalar const* values;

        PetscCallAbort(PETSC_COMM_WORLD,
                       MatGetRow(M_raw, i, &cols, &col_indices, &values));

        auto const columns = std::span<PetscInt const>(col_indices, cols);
        auto const values_M = std::span<PetscScalar const>(values, cols);

        for (auto const [j, mij] : std::views::zip(columns, values_M))
        {
            PetscScalar const fij = -mij * (xdot->get(j) - xdot->get(i));
            PetscCallAbort(PETSC_COMM_WORLD,
                           MatSetValue(F, i, j, fij, ADD_VALUES));
        }
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatRestoreRow(M_raw, i, &cols, &col_indices, &values));
    }
    finalize(F);
    auto P_plus = newZeroedInstance<GlobalVector>(matrix_specification);
    auto P_minus = newZeroedInstance<GlobalVector>(matrix_specification);

    P_plus->setLocalAccessibleVector();
    P_minus->setLocalAccessibleVector();
    // Fill P_plus, P_minus, Q_plus, Q_minus
    for (PetscInt i = start; i < end; ++i)
    {
        PetscInt cols;
        PetscInt const* col_indices;
        PetscScalar const* values;

        PetscCallAbort(PETSC_COMM_WORLD,
                       MatGetRow(F, i, &cols, &col_indices, &values));

        auto const column_indices =
            std::span<PetscInt const>(col_indices, cols);
        auto const values_span = std::span<PetscScalar const>(values, cols);
        for (auto const [j, fij] : std::views::zip(column_indices, values_span))
        {
            // skip diagonal
            if (i == j)
            {
                continue;
            }
            P_plus->add(i, std::max(0., fij));
            P_minus->add(i, std::min(0., fij));
        }
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatRestoreRow(F, i, &cols, &col_indices, &values));
    }
    P_plus->finalizeAssembly();
    P_plus->setLocalAccessibleVector();
    P_minus->finalizeAssembly();
    P_minus->setLocalAccessibleVector();
    auto Q_plus = newZeroedInstance<GlobalVector>(matrix_specification);
    auto Q_minus = newZeroedInstance<GlobalVector>(matrix_specification);
    Q_plus->setLocalAccessibleVector();
    Q_minus->setLocalAccessibleVector();

    for (PetscInt i = start; i < end; ++i)
    {
        PetscInt cols;
        PetscInt const* col_indices;

        PetscScalar Q_plus_i = 0.0;
        PetscScalar Q_minus_i = 0.0;

        PetscCallAbort(PETSC_COMM_WORLD,
                       MatGetRow(F, i, &cols, &col_indices, nullptr));
        auto const column_indices =
            std::span<PetscInt const>(col_indices, cols);
        for (auto const j : column_indices)
        {
            // skip diagonal
            if (i == j)
            {
                continue;
            }
            PetscScalar const x_prev_i = x_prev[process_id]->get(i);
            PetscScalar const x_prev_j = x_prev[process_id]->get(j);

            Q_plus_i = std::max({Q_plus_i, 0., x_prev_j - x_prev_i});
            Q_minus_i = std::min({Q_minus_i, 0., x_prev_j - x_prev_i});
        }
        Q_plus->set(i, Q_plus_i);
        Q_minus->set(i, Q_minus_i);
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatRestoreRow(F, i, &cols, &col_indices, nullptr));
    }
    Q_plus->finalizeAssembly();
    Q_plus->setLocalAccessibleVector();
    Q_minus->finalizeAssembly();
    Q_minus->setLocalAccessibleVector();

    auto R_plus = newZeroedInstance<GlobalVector>(matrix_specification);
    R_plus->setLocalAccessibleVector();

    auto row_sums_M = newZeroedInstance<GlobalVector>(matrix_specification);
    Vec M_L = row_sums_M->getRawVector();
    PetscCallAbort(PETSC_COMM_WORLD, MatGetRowSum(M_raw, M_L));
    row_sums_M->finalizeAssembly();
    row_sums_M->setLocalAccessibleVector();

    for (auto k = R_plus->getRangeBegin(); k < R_plus->getRangeEnd(); ++k)
    {
        PetscScalar const P_plus_i = P_plus->get(k);
        if (P_plus_i == 0.0)
        {
            continue;
        }
        PetscScalar const mi = row_sums_M->get(k);
        PetscScalar const Q_plus_i = Q_plus->get(k);
        R_plus->set(k, std::min(1.0, mi * Q_plus_i / dt / P_plus_i));
    }
    R_plus->finalizeAssembly();
    R_plus->setLocalAccessibleVector();
    auto R_minus = newZeroedInstance<GlobalVector>(matrix_specification);
    R_minus->setLocalAccessibleVector();

    for (auto k = R_minus->getRangeBegin(); k < R_minus->getRangeEnd(); ++k)
    {
        PetscScalar const P_minus_i = P_minus->get(k);
        if (P_minus_i == 0.0)
        {
            continue;
        }
        PetscScalar const mi = row_sums_M->get(k);
        PetscScalar const Q_minus_i = Q_minus->get(k);
        R_minus->set(k, std::min(1.0, mi * Q_minus_i / dt / P_minus_i));
    }
    R_minus->finalizeAssembly();
    R_minus->setLocalAccessibleVector();
    // walk over F, R_plus, and R_minus and compute alpha values; set entries in
    // the rhs b that limit the antidiffusive flux
    for (PetscInt i = start; i < end; ++i)
    {
        PetscInt cols;
        PetscInt const* col_indices;
        PetscScalar const* values;

        PetscCallAbort(PETSC_COMM_WORLD,
                       MatGetRow(F, i, &cols, &col_indices, &values));

        auto col_indices_span = std::span<PetscInt const>(col_indices, cols);
        auto values_span = std::span<PetscScalar const>(values, cols);

        double alpha_ij;

        for (auto [j, fij] : std::views::zip(col_indices_span, values_span))
        {
            if (i == j)
            {
                continue;
            }
            if (fij > 0.)
            {
                alpha_ij = std::min(R_plus->get(i), R_minus->get(j));
            }
            else
            {
                alpha_ij = std::min(R_minus->get(i), R_plus->get(j));
            }
            b.add(i, alpha_ij * fij);
        }
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatRestoreRow(F, i, &cols, &col_indices, &values));
    }

    // compute low-order operator
    // update K_raw with D: K_raw += D;
    for (PetscInt i = start; i < end; ++i)
    {
        PetscInt cols_D;
        PetscInt const* col_indices_D;
        PetscScalar const* values_D;
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatGetRow(D, i, &cols_D, &col_indices_D, &values_D));

        auto const column_indices =
            std::span<PetscInt const>(col_indices_D, cols_D);
        auto const values = std::span<PetscScalar const>(values_D, cols_D);
        for (auto [j, Dij] : std::views::zip(column_indices, values))
        {
            PetscCallAbort(PETSC_COMM_WORLD,
                           MatSetValues(K_raw, 1, &i, 1, &j, &Dij, ADD_VALUES));
        }
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatRestoreRow(D, i, &cols_D, &col_indices_D, &values_D));
    }
    finalize(K_raw);
    PetscCallAbort(PETSC_COMM_WORLD, MatScale(K_raw, -1.0));
    K.finalizeAssembly();
    M.setZero();
    for (PetscInt i = start; i < end; ++i)
    {
        PetscScalar const row_sum_M_i = row_sums_M->get(i);
        MatSetValues(M.getRawMatrix(), 1, &i, 1, &i, &row_sum_M_i,
                     INSERT_VALUES);
    }
    M.finalizeAssembly();
    MatDestroy(&D);
    MatDestroy(&F);
}
#endif  // end of ifdef USE_PETSC

void calculateFluxCorrectedTransport(
    const double t, const double dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    const MathLib::MatrixSpecifications& matrix_specification, GlobalMatrix& M,
    GlobalMatrix& K, GlobalVector& b)
{
#ifndef USE_PETSC
    calculateFluxCorrectedTransportSerial(t, dt, x, x_prev, process_id,
                                          matrix_specification, M, K, b);
#endif  // end of ifndef USE_PETSC

#ifdef USE_PETSC
    calculateFluxCorrectedTransportPETSc(t, dt, x, x_prev, process_id,
                                         matrix_specification, M, K, b);
#endif  // end of ifdef USE_PETSC
}

}  // namespace detail

void computeFluxCorrectedTransport(
    NumericalStabilization const& stabilizer, const double t, const double dt,
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    const MathLib::MatrixSpecifications& matrix_specification, GlobalMatrix& M,
    GlobalMatrix& K, GlobalVector& b)
{
    // if needed, calling the Flux-Corrected-Transport function
    return std::visit(
        [&](auto&& stabilizer)
        {
            using Stabilizer = std::decay_t<decltype(stabilizer)>;
            if constexpr (std::is_same_v<Stabilizer,
                                         NumLib::FluxCorrectedTransport>)
            {
                return detail::calculateFluxCorrectedTransport(
                    t, dt, x, x_prev, process_id, matrix_specification, M, K,
                    b);
            }
        },
        stabilizer);
}
}  // namespace NumLib
