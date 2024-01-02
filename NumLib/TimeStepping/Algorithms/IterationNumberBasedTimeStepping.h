/**
 * \file
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "TimeStepAlgorithm.h"

namespace NumLib
{
/**
 * \brief Iteration number based adaptive time stepping.
 *
 * This algorithm estimates a time step size depending on the number of
 * iterations (e.g. of iterative linear solvers, nonlinear methods, partitioned
 * coupling) needed in the last time step (see Hoffmann (2010) for
 * Newton-Raphson case).
 * The new time step \f$\Delta t_{n+1}\f$ size is calculated as
 * \f[
 *  \Delta t_{n+1} = \alpha \Delta t_n
 * \f]
 * with the previous time step size \f$\Delta t_{n}\f$ and a multiplier
 * coefficient
 * \f$\alpha\f$ depending on the iteration number.
 * Note that a time step size is always bounded by the minimum and maximum
 * allowed
 * value.
 * \f[
 *  \Delta t_{\min} \le \Delta t \le \Delta t_{\max}
 * \f]
 *
 * For example, users can setup the following time stepping strategy based on
 * the iteration number of the Newton-Raphson method in the last time step.
 * <table border="1">
 * <tr><th>Num. of Newton
 * steps</th><th>0-2</th><th>3-6</th><th>7-8</th><th>9<</th></tr>
 * <tr><th>Time step size
 * multiplier</th><th>1.6</th><th>1.</th><th>0.5</th><th>0.25 (repeat time
 * step)</th></tr>
 * <tr><th>Upper and lower bound</th><th colspan="4"> \f$ 1. \le \Delta t \le
 * 10.\f$ </th></tr>
 * </table>
 * A time step size is increased for the small iteration number, and decreased
 * for the
 * large iteration number. If the iteration number exceeds a user-defined
 * threshold (e.g. 9),
 * a time step is repeated with a smaller time step size.
 *
 * Reference
 * - Hoffmann J (2010) Reactive Transport and Mineral Dissolution/Precipitation
 * in Porous Media:Efficient Solution Algorithms, Benchmark Computations and
 * Existence of Global Solutions. PhD thesis. pp82.
 * Friedrich-Alexander-Universität Erlangen-Nürnberg.
 *
 */
class IterationNumberBasedTimeStepping final : public TimeStepAlgorithm
{
public:
    /**
     * @param t_initial             start time
     * @param t_end                 end time
     * @param min_dt                the minimum allowed time step size
     * @param max_dt                the maximum allowed time step size
     * @param initial_dt            initial time step size
     * @param iter_times_vector     a vector of iteration numbers
     * (\f$i_1\f$, \f$i_2\f$, ..., \f$i_n\f$) which defines intervals as
     * \f$[i_1,i_2)\f$, \f$[i_2,i_3)\f$, ..., \f$[i_n,\infty)\f$.
     * If an iteration number is larger than \f$i_n\f$, current time step is
     * repeated with the new time step size.
     * @param multiplier_vector     a vector of multiplier coefficients
     * (\f$a_1\f$, \f$a_2\f$, ..., \f$a_n\f$) corresponding to the intervals
     * given by iter_times_vector.
     * A time step size is calculated by \f$\Delta t_{n+1} = a * \Delta t_{n}\f$
     * @param fixed_times_for_output a vector of fixed time points for output
     */
    IterationNumberBasedTimeStepping(
        double const t_initial,
        double const t_end,
        double const min_dt,
        double const max_dt,
        double const initial_dt,
        std::vector<int>&& iter_times_vector,
        std::vector<double>&& multiplier_vector,
        std::vector<double> const& fixed_times_for_output);

    ~IterationNumberBasedTimeStepping() override = default;

    std::tuple<bool, double> next(double solution_error, int number_iterations,
                                  NumLib::TimeStep& ts_previous,
                                  NumLib::TimeStep& ts_current) override;

    bool isSolutionErrorComputationNeeded() const override { return true; }

    bool canReduceTimestepSize(
        NumLib::TimeStep const& timestep_previous,
        NumLib::TimeStep const& timestep_current) const override;

private:
    /// Calculate the next time step size.
    double getNextTimeStepSize(NumLib::TimeStep const& ts_previous,
                               NumLib::TimeStep const& ts_current) const;

    /// Find a multiplier for the given number of iterations.
    double findMultiplier(int const number_iterations,
                          NumLib::TimeStep const& ts_current) const;

    /// This vector stores the number of iterations to which the respective
    /// multiplier coefficient will be applied.
    const std::vector<int> _iter_times_vector;
    /// This vector stores the multiplier coefficients.
    const std::vector<double> _multiplier_vector;
    /// The minimum allowed time step size.
    const double _min_dt;
    /// The maximum allowed time step size.
    const double _max_dt;
    /// Initial time step size.
    const double _initial_dt;
    /// The maximum allowed iteration number to accept current time step.
    const int _max_iter;
    /// The number of nonlinear iterations.
    int _iter_times = 0;

    bool _previous_time_step_accepted = true;
    std::vector<double> const _fixed_times_for_output;
};

}  // namespace NumLib
