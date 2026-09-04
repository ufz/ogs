// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "TimeStepAlgorithm.h"

namespace NumLib
{
using RepeatDtPair = std::tuple<std::size_t, double>;

/**
 * \brief Fixed time stepping algorithm
 *
 * This algorithm returns time step size defined by a user priori.
 *
 * If the prescribed time step sizes are used up before the end time is
 * reached---which happens if the steps actually taken are smaller than the
 * prescribed ones, as in the staggered coupling scheme---the last prescribed
 * time step size taken by the scheme is repeated until the end time.
 */
class FixedTimeStepping final : public TimeStepAlgorithm
{
public:
    /**
     * Constructor with homogeneous time step size
     *
     * A user provides a single time step size \f$\Delta t\f$. Total number of
     * time steps is calculated by
     * \f[
     *  n=\frac{t_{\rm n} - t_0}{\Delta t}
     * \f].
     *
     * \pre \c t_end is larger than \c t0 and \c dt is positive; both are
     * required for a non-empty sequence of time steps and are checked with
     * OGS_FATAL.
     *
     * @param t0 start time
     * @param t_end end time
     * @param dt uniform time step size
     */
    FixedTimeStepping(double t0, double t_end, double dt);

    /**
     * Constructor with user-specified time step sizes including additional time
     * steps for output.
     *
     * \pre \c tn is larger than \c t0 and \c repeat_dt_pairs is valid in the
     * sense of areRepeatDtPairsValid(), i.e. non-empty with positive
     * \c delta_t and non-zero \c repeat entries; both are required for a
     * non-empty sequence of time steps and are checked with OGS_FATAL.
     *
     * @param t0 start time
     * @param tn end time
     * @param repeat_dt_pairs pairs of the number of repetitions and the time
     * step size
     * @param fixed_times_for_output times at which output is requested; the
     * time step containing such a time is split at it
     */
    FixedTimeStepping(double t0, double tn,
                      std::vector<RepeatDtPair> const& repeat_dt_pairs,
                      std::vector<double> const& fixed_times_for_output);

    double next(double solution_error, int number_iterations,
                NumLib::TimeStep& ts_previous,
                NumLib::TimeStep& ts_current) override;

    static bool areRepeatDtPairsValid(
        std::vector<RepeatDtPair> const& repeat_dt_pairs);

private:
    /// a vector of time step sizes
    std::vector<double> dt_vector_;

    /// The size of the last time step the scheme takes as prescribed by the
    /// user. Kept separately from dt_vector_ because the latter's last entry is
    /// not necessarily a prescribed size: incorporateFixedTimesForOutput()
    /// splits the interval containing a fixed output time into two smaller
    /// ones. It is not the last (repeat, delta_t) pair's size either, because
    /// pairs whose steps start beyond the end time are not part of the scheme.
    double last_prescribed_dt_ = 0.0;

    /// set if it has already been reported that the prescribed time step sizes
    /// are used up and the last one is repeated, in order to report it once
    /// only.
    bool reuse_of_last_dt_reported_ = false;
};

std::size_t findDeltatInterval(Time const& t_initial,
                               std::vector<double> const& delta_ts,
                               Time const& fixed_output_time);
}  // namespace NumLib
