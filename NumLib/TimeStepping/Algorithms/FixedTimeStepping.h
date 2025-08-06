/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
     * @param t0 start time
     * @param t_end end time
     * @param dt uniform time step size
     */
    FixedTimeStepping(double t0, double t_end, double dt);

    /**
     * Constructor with user-specified time step sizes including additional time
     * steps for output.
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
    std::vector<double> _dt_vector;
};

std::size_t findDeltatInterval(Time const& t_initial,
                               std::vector<double> const& delta_ts,
                               Time const& fixed_output_time);
}  // namespace NumLib
