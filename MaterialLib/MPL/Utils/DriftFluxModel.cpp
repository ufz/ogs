// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DriftFluxModel.h"

#include <Eigen/Dense>
#include <cmath>

#include "BaseLib/Logging.h"
#include "NumLib/NewtonRaphson.h"

namespace MaterialPropertyLib
{
static double voidFractionResidual(double const alpha, double const dryness,
                                   double const vapour_water_density,
                                   double const liquid_water_density,
                                   double const v_mix, double const C_0,
                                   double const u_gu)
{
    double const rho_mix =
        alpha * vapour_water_density + (1 - alpha) * liquid_water_density;

    return dryness * liquid_water_density * rho_mix * v_mix -
           alpha * C_0 * dryness * liquid_water_density * rho_mix * v_mix -
           alpha * C_0 * (1 - dryness) * vapour_water_density * rho_mix *
               v_mix -
           alpha * vapour_water_density * liquid_water_density * u_gu;
}

static double voidFractionResidualDerivative(
    double const alpha, double const dryness, double const vapour_water_density,
    double const liquid_water_density, double const v_mix, double const C_0,
    double const u_gu)
{
    return dryness * liquid_water_density * v_mix *
               (vapour_water_density - liquid_water_density) -
           (C_0 * dryness * liquid_water_density +
            C_0 * (1 - dryness) * vapour_water_density) *
               (2 * alpha * vapour_water_density +
                (1 - 2 * alpha) * liquid_water_density) *
               v_mix -
           vapour_water_density * liquid_water_density * u_gu;
}

double driftFluxProfileParameter(double const dryness)
{
    return 1 + 0.12 * (1 - dryness);
}

double driftFluxVelocity(double const dryness, double const temperature,
                         double const vapour_water_density,
                         double const liquid_water_density)
{
    double const sigma_gl = 0.2358 *
                            std::pow((1 - temperature / 647.096), 1.256) *
                            (1 - 0.625 * (1 - temperature / 647.096));

    return 1.18 * (1 - dryness) *
           std::pow((9.81) * sigma_gl *
                        (liquid_water_density - vapour_water_density),
                    0.25) /
           std::pow(liquid_water_density, 0.5);
}

double computeVapourVoidFraction(double const dryness,
                                 double const vapour_water_density,
                                 double const liquid_water_density,
                                 double const v_mix, double const C_0,
                                 double const u_gu)
{
    double alpha = 0;

    if (dryness == 0)
    {
        return alpha;
    }

    using LocalJacobianMatrix = Eigen::Matrix<double, 1, 1, Eigen::RowMajor>;
    using LocalResidualVector = Eigen::Matrix<double, 1, 1>;
    using LocalUnknownVector = Eigen::Matrix<double, 1, 1>;
    LocalJacobianMatrix J_loc;

    Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(1);

    auto const update_residual = [&](LocalResidualVector& residual)
    {
        residual(0) =
            voidFractionResidual(alpha, dryness, vapour_water_density,
                                 liquid_water_density, v_mix, C_0, u_gu);
    };

    auto const update_jacobian = [&](LocalJacobianMatrix& jacobian)
    {
        jacobian(0) = voidFractionResidualDerivative(
            alpha, dryness, vapour_water_density, liquid_water_density, v_mix,
            C_0, u_gu);
    };

    auto const update_solution = [&](LocalUnknownVector const& increment)
    { alpha += increment[0]; };

    const int maximum_iterations(20);
    const double residuum_tolerance(1.e-10);
    const double increment_tolerance(0);

    auto newton_solver = NumLib::NewtonRaphson(
        linear_solver, update_jacobian, update_residual, update_solution,
        {maximum_iterations, residuum_tolerance, increment_tolerance});

    auto const success_iterations = newton_solver.solve(J_loc);

    if (!success_iterations)
    {
        WARN(
            "Attention! Steam void fraction has not been correctly "
            "calculated!");
    }

    return alpha;
}

double mixtureSlipParameter(double const alpha,
                            double const vapour_water_density,
                            double const liquid_water_density,
                            double const v_mix, double const C_0,
                            double const u_gu)
{
    // The zero slip limit at a void fraction of one, see the documentation of
    // this function.
    if (alpha == 1)
    {
        return 0.;
    }

    double const rho_mix =
        alpha * vapour_water_density + (1 - alpha) * liquid_water_density;

    return alpha * liquid_water_density * vapour_water_density * rho_mix /
           (1 - alpha) /
           std::pow((alpha * C_0 * vapour_water_density +
                     (1 - alpha * C_0) * liquid_water_density),
                    2) *
           std::pow((C_0 - 1) * v_mix + u_gu, 2);
}
}  // namespace MaterialPropertyLib
