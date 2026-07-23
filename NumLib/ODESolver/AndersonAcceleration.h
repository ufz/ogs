// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

namespace NumLib::detail
{
/*! Computes the Anderson mixing weights \f$ \theta \f$ minimising
 * \f$ \|\sum_i \theta_i f_i\|^2 \f$ subject to \f$ \sum_i \theta_i = 1 \f$,
 * given the Gram matrix \f$ G = F^T F \f$ of the stored residuals
 * \f$ f_i = g(x_i) - x_i \f$.
 *
 * The constrained least-squares problem is solved via the Lagrange
 * formulation
 * \f[
 *   \begin{bmatrix} G & \mathbf{1} \\ \mathbf{1}^T & 0 \end{bmatrix}
 *   \begin{bmatrix} \theta \\ \lambda \end{bmatrix}
 *   =
 *   \begin{bmatrix} \mathbf{0} \\ 1 \end{bmatrix}
 * \f]
 *
 * \p G is normalized internally by its largest diagonal entry
 * (\f$ \max_i \|f_i\|^2 \f$) to keep the augmented system well-conditioned
 * regardless of the residuals' magnitude; the solution \f$ \theta \f$ is
 * invariant under this scaling, because the scale factor cancels in
 * \f$ \theta = G^{-1}\mathbf{1} / (\mathbf{1}^T G^{-1}\mathbf{1}) \f$.
 *
 * Linearly dependent (or nearly dependent) residuals are the classical failure
 * mode of Anderson acceleration. The computed mixture is therefore accepted
 * only if
 * -# it predicts a smaller residual norm than the plain step it would replace,
 *    which rules out the arbitrary minimizer that a degenerate history admits,
 *    and
 * -# its weights stay bounded, which rules out the mixtures of near-identical
 *    iterates whose value is pure cancellation error.
 *
 * Otherwise \f$ \theta = (0,\dots,0,1) \f$ is returned, i.e. unit weight on the
 * newest stored step, which reproduces the plain (damped) Picard update. The
 * same fallback applies when all stored residuals vanish. Acceleration thus
 * degrades to plain Picard instead of amplifying rounding error.
 *
 * \pre \p G is symmetric positive semi-definite and at least 1x1.
 */
Eigen::VectorXd computeAndersonWeights(Eigen::MatrixXd G);

}  // namespace NumLib::detail
