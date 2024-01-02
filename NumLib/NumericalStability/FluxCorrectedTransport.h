/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MathLib/LinAlg/MatrixSpecifications.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumericalStabilization.h"

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

void calculateFluxCorrectedTransport(
    [[maybe_unused]] const double t, [[maybe_unused]] const double dt,
    [[maybe_unused]] std::vector<GlobalVector*> const& x,
    [[maybe_unused]] std::vector<GlobalVector*> const& x_prev,
    [[maybe_unused]] int const process_id,
    [[maybe_unused]] const MathLib::MatrixSpecifications& matrix_specification,
    [[maybe_unused]] GlobalMatrix& M, [[maybe_unused]] GlobalMatrix& K,
    [[maybe_unused]] GlobalVector& b)
{
#ifndef USE_PETSC
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
        double const mi = M_L(k);
        double const Q_plus_i = Q_plus->get(k);
        double const P_plus_i = P_plus->get(k);

        double const tmp =
            P_plus_i == 0. ? 0.0 : std::min(1.0, mi * Q_plus_i / dt / P_plus_i);

        R_plus->set(k, tmp);
    }

    for (auto k = R_minus->getRangeBegin(); k < R_minus->getRangeEnd(); ++k)
    {
        double const mi = M_L(k);
        double const Q_minus_i = Q_minus->get(k);
        double const P_minus_i = P_minus->get(k);

        double const tmp = P_minus_i == 0.
                               ? 0.0
                               : std::min(1.0, mi * Q_minus_i / dt / P_minus_i);
        R_minus->set(k, tmp);
    }

    auto alpha = newZeroedInstance<GlobalMatrix>(matrix_specification);
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

                alpha->setValue(it.row(), it.col(), alpha_ij);
            }
            else
            {
                double const R_minus_i = R_minus->get(it.row());
                double const R_plus_j = R_plus->get(it.col());
                double const alpha_ij = std::min(R_minus_i, R_plus_j);

                alpha->setValue(it.row(), it.col(), alpha_ij);
            }
        }
    }

    // compute limited antidiffusive fluxes
    for (int k = 0; k < F_raw.outerSize(); ++k)
    {
        for (RawMatrixType::InnerIterator it(F_raw, k); it; ++it)
        {
            if (it.row() != it.col())
            {
                double const fij = it.value();
                double const alpha_ij = alpha->get(it.row(), it.col());

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
#endif  // end of ifndef USE_PETSC
}
}  // namespace detail

inline void computeFluxCorrectedTransport(
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
