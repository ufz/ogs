// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AndersonAcceleration.h"

#include <spdlog/fmt/ranges.h>

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>

#include "AndersonWeights.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"

namespace NumLib
{
namespace
{
//! Whether the differences f_i - f_newest of the steps with Gram matrix \p G
//! are linearly independent, i.e. whether the constrained least-squares
//! problem has a unique minimiser. A difference shorter than sqrt(eps) times
//! the largest step counts as a duplicate; otherwise the differences count as
//! dependent if the smallest eigenvalue of their unit-diagonal Gram matrix
//! (their angles) is below sqrt(eps).
bool hasUniqueMinimizer(Eigen::MatrixXd const& G)
{
    int const n = static_cast<int>(G.rows()) - 1;  // number of differences
    if (n < 1)
    {
        return true;
    }
    double const tolerance = std::sqrt(std::numeric_limits<double>::epsilon());

    // minCoeff() skips a NaN unless it is the first entry, hence hasNaN().
    auto const allAboveTolerance =
        [](Eigen::VectorXd const& v, double const threshold)
    { return !v.hasNaN() && v.minCoeff() > threshold; };

    // Gram matrix of the differences d_i = f_i - f_n.
    Eigen::MatrixXd D = G.topLeftCorner(n, n);
    D.colwise() -= G.col(n).head(n);
    D.rowwise() -= G.row(n).head(n);
    D.array() += G(n, n);

    // A difference negligible against the largest step is a duplicate step.
    Eigen::VectorXd const norms = D.diagonal().cwiseMax(0.0).cwiseSqrt();
    if (!allAboveTolerance(norms,
                           tolerance * std::sqrt(G.diagonal().maxCoeff())))
    {
        return false;
    }
    Eigen::VectorXd const inv_norms = norms.cwiseInverse();
    Eigen::MatrixXd const C =
        inv_norms.asDiagonal() * D * inv_norms.asDiagonal();
    Eigen::VectorXd const eigenvalues =
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(C,
                                                       Eigen::EigenvaluesOnly)
            .eigenvalues();
    return allAboveTolerance(eigenvalues, tolerance);
}
}  // namespace

namespace detail
{
Eigen::VectorXd computeAndersonWeights(Eigen::MatrixXd G)
{
    int const history_size = static_cast<int>(G.rows());

    // At least one residual is required: G is the Gram matrix of the stored
    // steps, so maxCoeff() below is a reduction over a non-empty diagonal.
    if (history_size < 1)
    {
        OGS_FATAL(
            "Anderson acceleration: the mixing weights were requested for an "
            "empty history (Gram matrix of size {:d}x{:d}).",
            history_size, static_cast<int>(G.cols()));
    }

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

    // Drop the oldest steps until the minimiser is unique (e.g. duplicate
    // steps, or a map acting on one DoF).
    int number_of_steps = history_size;
    while (number_of_steps > 1 && !hasUniqueMinimizer(G.bottomRightCorner(
                                      number_of_steps, number_of_steps)))
    {
        --number_of_steps;
    }
    if (number_of_steps == 1)
    {
        // Only the plain step is left; the descent test would wrongly
        // reject it.
        return newest_step_only();
    }
    // Renormalised: the dropped steps may have been much larger.
    Eigen::MatrixXd G_used =
        G.bottomRightCorner(number_of_steps, number_of_steps);
    G_used /= G_used.diagonal().maxCoeff();

    Eigen::MatrixXd M(number_of_steps + 1, number_of_steps + 1);
    M.topLeftCorner(number_of_steps, number_of_steps) = G_used;
    M.topRightCorner(number_of_steps, 1).setOnes();
    M.bottomLeftCorner(1, number_of_steps).setOnes();
    M(number_of_steps, number_of_steps) = 0.0;

    Eigen::VectorXd rhs_aa(number_of_steps + 1);
    rhs_aa.setZero();
    rhs_aa(number_of_steps) = 1.0;

    Eigen::VectorXd theta = Eigen::VectorXd::Zero(history_size);
    theta.tail(number_of_steps) =
        M.fullPivLu().solve(rhs_aa).head(number_of_steps);

    if (!theta.allFinite())
    {
        // A fallback silently changes the iterate the solver would otherwise
        // take, so it is reported at INFO level rather than hidden in DBUG.
        INFO(
            "Anderson acceleration: the mixing weights came out non-finite. "
            "Falling back to the plain Picard step for this iteration.");
        return newest_step_only();
    }

    // Unique weights may still be useless; validate them.
    //
    // (a) Descent: the mixture only earns its place if the residual norm it
    //     predicts is smaller than that of the plain step it would replace.
    //     Both norms are read off the same normalized G (divided by g_scale
    //     above), so the shared 1/g_scale factor cancels and the comparison is
    //     exactly the one on the unscaled residuals.
    double const mixed_residual_norm_2 = theta.dot(G * theta);
    double const newest_step_norm_2 = G(history_size - 1, history_size - 1);

    // (b) No long lever arms: weights far outside [0, 1] mean the mixed
    //     iterate is a difference of near-identical vectors, i.e. dominated by
    //     cancellation. Healthy histories stay at |theta| ~ 1, whereas steps
    //     that agree to k digits produce weights of magnitude 10^k.
    //
    //     The threshold is a heuristic, not a derived bound. A weight of
    //     magnitude 10^k sacrifices about k of the ~16 significant decimal
    //     digits of a double to cancellation; capping at 10^2 admits the
    //     modest lever arms of genuinely useful mixing (empirically |theta| up
    //     to ~10) while rejecting the 10^3-and-up weights that signal a
    //     degenerate history. It is deliberately loose: the descent test (a)
    //     is the primary guard, and this one only catches the cancellation
    //     cases that slip past it.
    constexpr double max_weight = 1e2;

    // (c) A numerically zero predicted residual is exact only for a map that
    //     is affine on the span of the steps; otherwise it is a secant
    //     extrapolation of length ~|theta| steps (f_1 = 1.02 f_0 gives
    //     theta = (51, -50)). Accept it only for |theta| <= 10. "Zero" is
    //     relative to the rounding error of theta^T G theta, which grows like
    //     eps ||theta||_1^2 max G_ii; sqrt(eps) gives a wide margin.
    constexpr double max_exact_extrapolation_weight = 1e1;
    double const theta_l1 = theta.lpNorm<1>();
    double const used_scale = G.diagonal().tail(number_of_steps).maxCoeff();
    double const numerically_zero_2 =
        std::sqrt(std::numeric_limits<double>::epsilon()) * theta_l1 *
        theta_l1 * used_scale;
    bool const long_exact_extrapolation =
        mixed_residual_norm_2 <= numerically_zero_2 &&
        theta.cwiseAbs().maxCoeff() > max_exact_extrapolation_weight;

    // Convergence is checked on the unmixed step (see NonlinearSolver), so a
    // rejected mixture costs an iteration, not correctness.
    //
    // Negated comparison, so that a NaN norm takes the fallback as well.
    if (!(mixed_residual_norm_2 < newest_step_norm_2) ||
        theta.cwiseAbs().maxCoeff() > max_weight || long_exact_extrapolation)
    {
        // A fallback silently changes the iterate the solver would otherwise
        // take, so it is reported at INFO level rather than hidden in DBUG.
        INFO(
            "Anderson acceleration: rejected the mixture of {:d} stored steps "
            "(predicted residual {:g} vs. {:g} for the plain step, largest "
            "weight {:g}). Falling back to the plain Picard step for this "
            "iteration.",
            number_of_steps, mixed_residual_norm_2, newest_step_norm_2,
            theta.cwiseAbs().maxCoeff());
        return newest_step_only();
    }

    return theta;
}

}  // namespace detail

AndersonAcceleration::AndersonAcceleration(int const depth)
    : _depth(depth), _gram(depth, depth)
{
    if (_depth >= min_mixing_depth)
    {
        _history.reserve(_depth);
    }
}

AndersonAcceleration::~AndersonAcceleration()
{
    for (auto const& entry : _history)
    {
        releaseHistoryEntry(entry);
    }
}

void AndersonAcceleration::releaseHistoryEntry(HistoryEntry const& entry)
{
    NumLib::GlobalVectorProvider::provider.releaseVector(*entry.x);
    NumLib::GlobalVectorProvider::provider.releaseVector(*entry.f);
}

void AndersonAcceleration::accelerate(GlobalVector const& x_old,
                                      GlobalVector& x_new)
{
    namespace LinAlg = MathLib::LinAlg;

    // A depth below min_mixing_depth admits no mixing (plain Picard); nothing
    // is stored.
    if (_depth < min_mixing_depth)
    {
        return;
    }

    // Additionally mixes the last _depth damped steps
    //   f_i = beta*(g(x_i) - x_i) (i.e. x_new - x_old computed after the beta
    //   relaxation already applied by the caller) to find the optimal theta
    //   minimising ||sum theta_i f_i|| s.t. sum theta_i = 1, then sets
    //   x_new = sum theta_i*(x_i + f_i).

    // Whether the circular buffer is full and the oldest entry is about to be
    // evicted (needed for the incremental Gram update).
    bool const rotated = static_cast<int>(_history.size()) == _depth;
    if (!rotated)
    {
        // The id out-params are unused: the provider allocates a fresh vector
        // on every call and never re-fetches by id.
        std::size_t x_id = 0u;
        std::size_t f_id = 0u;
        _history.push_back(
            {&NumLib::GlobalVectorProvider::provider.getVector(x_id),
             &NumLib::GlobalVectorProvider::provider.getVector(f_id)});
    }
    else
    {
        // Recycle the oldest entry as the newest one.
        std::rotate(_history.begin(), _history.begin() + 1, _history.end());
    }

    auto const& newest = _history.back();

    // x = x_old, f = x_new - x_old
    LinAlg::copy(x_old, *newest.x);
    LinAlg::copy(x_new, *newest.f);
    LinAlg::axpy(*newest.f, -1.0, x_old);

    // Actual window size, <= _depth while the buffer fills.
    int const history_size = static_cast<int>(_history.size());

    // Incrementally maintain the (history_size x history_size) Gram matrix
    // G = F^T F whose columns are the stored damped steps f_0 ...
    // f_{history_size-1}. All steps but the newest are unchanged from the
    // previous iteration, so only the last row/column is recomputed -
    // history_size dot products instead of a full
    // history_size*(history_size+1)/2 rebuild. On a rotate the oldest entry
    // (index 0) was evicted, so the cached block is first shifted up-left by
    // one.
    if (rotated)
    {
        _gram.topLeftCorner(history_size - 1, history_size - 1) =
            _gram.block(1, 1, history_size - 1, history_size - 1).eval();
    }
    for (int i = 0; i < history_size; ++i)
    {
        double const d = LinAlg::dot(*_history[i].f, *newest.f);
        _gram(i, history_size - 1) = d;
        _gram(history_size - 1, i) = d;
    }

    // A single stored step needs no mixing: the sum-to-one constraint forces
    // theta = (1), which just reproduces the damped step already held in x_new.
    if (history_size < min_mixing_depth)
    {
        return;
    }

    // Solve G theta = e (least-squares) with the constraint sum theta_i = 1 via
    // a simple Lagrange formulation:
    //
    //   [ G  1 ] [ theta  ] = [ 0 ]
    //   [ 1  0 ] [ lambda ]   [ 1 ]
    //
    // The beta factor scales G by beta^2 and cancels in theta, so the weights
    // are identical to the undamped case. The Anderson update is then:
    //   x_anderson = sum_i theta_i * (x_i + f_i)
    //              = sum_i theta_i * (x_i + beta*(g(x_i)-x_i))
    //              = sum_i theta_i * ((1-beta)*x_i + beta*g(x_i))
    Eigen::MatrixXd const G = _gram.topLeftCorner(history_size, history_size);

    Eigen::VectorXd const theta = detail::computeAndersonWeights(G);

    // Accumulate the Anderson mixed iterate directly into x_new. Its previous
    // value is no longer needed: the newest step was already extracted from it
    // above, and it is not aliased by any history entry (those are independent
    // copies).
    x_new.setZero();
    for (int i = 0; i < history_size; ++i)
    {
        // x_new += theta_i * (x_i + f_i)
        LinAlg::axpy(x_new, theta(i), *_history[i].x);
        LinAlg::axpy(x_new, theta(i), *_history[i].f);
    }

    DBUG("Picard/Anderson: history size {:d}, theta=[{:.4g}]", history_size,
         fmt::join(theta.data(), theta.data() + history_size, ", "));
}

void AndersonAcceleration::dropLastStep()
{
    // Drop the just-added (now stale) history entry, since the iteration is
    // being repeated. In the full-buffer case this is the recycled slot;
    // releasing it by reference is safe because the provider tracks vectors by
    // pointer, not by the (unused) id.
    //
    // The rotation and the Gram shift performed by accelerate() are not undone,
    // and need not be: dropping the newest entry leaves the buffer holding the
    // remaining entries in order, and the shifted top-left block of the Gram
    // matrix is exactly their Gram matrix. The oldest entry stays evicted,
    // which merely shortens the sliding window by one.
    if (!_history.empty())
    {
        releaseHistoryEntry(_history.back());
        _history.pop_back();
    }
}

}  // namespace NumLib
