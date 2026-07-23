// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AndersonAcceleration.h"

#include <Eigen/LU>
#include <cassert>

#include "BaseLib/Logging.h"

namespace NumLib::detail
{
Eigen::VectorXd computeAndersonWeights(Eigen::MatrixXd G)
{
    int const history_size = static_cast<int>(G.rows());

    // At least one residual is required: G is the Gram matrix of the stored
    // steps, so maxCoeff() below is a reduction over a non-empty diagonal.
    assert(history_size >= 1);

    // Unit weight on the newest stored step, i.e. the plain (damped) Picard
    // update. Used whenever the stored history does not admit trustworthy
    // mixing weights.
    auto const newest_step_only = [history_size]()
    {
        Eigen::VectorXd theta = Eigen::VectorXd::Zero(history_size);
        theta(history_size - 1) = 1.0;
        return theta;
    };

    double const g_scale = G.diagonal().maxCoeff();
    // Negated comparison, so that a NaN scale takes this branch as well.
    if (!(g_scale > 0.0))
    {
        // All stored steps vanish; there is nothing left to mix.
        return newest_step_only();
    }
    G /= g_scale;

    Eigen::MatrixXd M(history_size + 1, history_size + 1);
    M.topLeftCorner(history_size, history_size) = G;
    M.topRightCorner(history_size, 1).setOnes();
    M.bottomLeftCorner(1, history_size).setOnes();
    M(history_size, history_size) = 0.0;

    Eigen::VectorXd rhs_aa(history_size + 1);
    rhs_aa.setZero();
    rhs_aa(history_size) = 1.0;

    Eigen::VectorXd const theta =
        M.fullPivLu().solve(rhs_aa).head(history_size);

    if (!theta.allFinite())
    {
        DBUG(
            "Anderson acceleration: the mixing weights came out non-finite. "
            "Falling back to the plain Picard step for this iteration.");
        return newest_step_only();
    }

    // The rank-revealing solve above never fails outright, so the weights have
    // to be validated on their own merits. Note that an invertibility test on
    // M would be the wrong check: it rejects healthy histories whose residual
    // norms span many orders of magnitude (the normal situation for a
    // converging iteration) while passing the ill-conditioned cases that
    // actually do harm.
    //
    // (a) Descent: the mixture only earns its place if the residual norm it
    //     predicts is smaller than that of the plain step it would replace.
    //     For a degenerate history - linearly dependent steps, in particular
    //     duplicates - the minimizer is not unique and the solve returns an
    //     arbitrary one, which this test discards.
    double const mixed_residual_norm_2 = theta.dot(G * theta);
    double const newest_step_norm_2 = G(history_size - 1, history_size - 1);

    // (b) No long lever arms: weights far outside [0, 1] mean the mixed
    //     iterate is a difference of near-identical vectors, i.e. dominated by
    //     cancellation. Healthy histories stay at |theta| ~ 1, whereas steps
    //     that agree to k digits produce weights of magnitude 10^k.
    constexpr double max_weight = 1e2;

    // Negated comparison, so that a NaN norm takes the fallback as well.
    if (!(mixed_residual_norm_2 < newest_step_norm_2) ||
        theta.cwiseAbs().maxCoeff() > max_weight)
    {
        DBUG(
            "Anderson acceleration: rejected the mixture of {:d} stored steps "
            "(predicted residual {:g} vs. {:g} for the plain step, largest "
            "weight {:g}). Falling back to the plain Picard step for this "
            "iteration.",
            history_size, mixed_residual_norm_2, newest_step_norm_2,
            theta.cwiseAbs().maxCoeff());
        return newest_step_only();
    }

    return theta;
}

}  // namespace NumLib::detail
