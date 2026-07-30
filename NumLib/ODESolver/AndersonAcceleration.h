// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace NumLib
{
/*! Anderson acceleration of the damped Picard fixpoint iteration.
 *
 * Owns a sliding window of the last \c depth damped steps and the Gram matrix
 * of that window, and mixes them into an accelerated iterate. History vectors
 * are taken from the global vector provider on demand and returned in the
 * destructor (RAII), so no manual release is required at the call site.
 *
 * The window and its Gram matrix are only ever mutated together; keeping them
 * in one object makes "history and Gram stay in sync" a class invariant instead
 * of a convention spread across the solver loop.
 *
 * A \c depth below \c min_mixing_depth admits no mixing (the sum-to-one
 * constraint forces unit weight on the single stored step), so such an instance
 * is an inert no-op: it takes no vector from the global vector provider and
 * leaves the iterate untouched.
 * This is the plain-Picard case.
 */
class AndersonAcceleration final
{
public:
    //! Smallest \c depth that admits mixing: below it the sum-to-one constraint
    //! forces unit weight on the single stored step, i.e. plain Picard.
    static constexpr int min_mixing_depth = 2;

    //! \param depth number of previous iterates retained for mixing; a value
    //!              below \c min_mixing_depth makes the instance an inert no-op
    //!              (plain Picard).
    explicit AndersonAcceleration(int depth);

    ~AndersonAcceleration();

    AndersonAcceleration(AndersonAcceleration const&) = delete;
    AndersonAcceleration& operator=(AndersonAcceleration const&) = delete;

    /*! Records the newest damped step \f$ x_{\rm old} \to x_{\rm new} \f$ and
     * overwrites \p x_new in place with the Anderson-mixed iterate.
     *
     * While fewer than two steps are stored, or when the mixture is rejected as
     * untrustworthy (see \c detail::computeAndersonWeights), \p x_new is left
     * unchanged. For a no-op instance (\c depth < \c min_mixing_depth) this
     * does nothing.
     *
     * \param x_old the iterate entering this step.
     * \param x_new in: the damped Picard step \f$ f = \beta(g(x_{\rm old}) -
     *              x_{\rm old}) \f$ added to \p x_old; out: the mixed iterate.
     */
    void accelerate(GlobalVector const& x_old, GlobalVector& x_new);

    //! Discards the step recorded by the most recent \c accelerate() call, used
    //! when the current iteration is repeated. The Gram shift performed while
    //! recording is deliberately not undone (see implementation).
    void dropLastStep();

private:
    //! One entry of the history: the iterate \c x it was taken from and the
    //! (possibly damped) step \f$ f = \beta(g(x) - x) \f$ leading away from it.
    //! The two vectors are only ever appended, rotated and dropped together.
    struct HistoryEntry
    {
        GlobalVector* x;
        GlobalVector* f;
    };

    //! Returns \p entry's vectors to the global vector provider.
    static void releaseHistoryEntry(HistoryEntry const& entry);

    //! Maximum window size; mixing is active only for values >=
    //! \c min_mixing_depth.
    int const _depth;

    //! Circular buffer of history entries, oldest first (size <= _depth).
    std::vector<HistoryEntry> _history;

    //! Gram matrix G = F^T F of the stored steps, maintained incrementally
    //! across iterations (only the newest step's row/column is recomputed).
    //! Sized once to the maximum window (_depth x _depth) so the incremental
    //! update never reallocates; only the leading history_size x history_size
    //! block is live while the window fills.
    Eigen::MatrixXd _gram;
};

}  // namespace NumLib
