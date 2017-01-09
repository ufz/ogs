/**
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "ITimeStepAlgorithm.h"

namespace NumLib
{

/**
 * \brief Iteration number based adaptive time stepping
 *
 * This algorithm estimates a time step size depending on the number of
 * iterations (e.g. of iterative linear solvers, nonlinear methods, partitioned
 * coupling) needed in the last time step (see Hoffmann (2010) for Newton-Raphson case).
 * The new time step \f$\Delta t_{n+1}\f$ size is calculated as
 * \f[
 *  \Delta t_{n+1} = \alpha \Delta t_n
 * \f]
 * with the previous time step size \f$\Delta t_{n}\f$ and a multiplier coefficient
 * \f$\alpha\f$ depending on the iteration number.
 * Note that a time step size is always bounded by the minimum and maximum allowed
 * value.
 * \f[
 *  \Delta t_{\min} \le \Delta t \le \Delta t_{\max}
 * \f]
 *
 * For example, users can setup the following time stepping strategy based on
 * the iteration number of the Newton-Raphson method in the last time step.
 * <table border="1">
 * <tr><th>Num. of Newton steps</th><th>0-2</th><th>3-6</th><th>7-8</th><th>9<</th></tr>
 * <tr><th>Time step size multiplier</th><th>1.6</th><th>1.</th><th>0.5</th><th>0.25 (repeat time step)</th></tr>
 * <tr><th>Upper and lower bound</th><th colspan="4"> \f$ 1. \le \Delta t \le 10.\f$ </th></tr>
 * </table>
 * A time step size is increased for the small iteration number, and decreased for the
 * large iteration number. If the iteration number exceeds a user-defined threshold (e.g. 9),
 * a time step is repeated with a smaller time step size.
 *
 * Reference
 * - Hoffmann J (2010) Reactive Transport and Mineral Dissolution/Precipitation
 * in Porous Media:Efficient Solution Algorithms, Benchmark Computations and
 * Existence of Global Solutions. PhD thesis. pp82. Friedrich-Alexander-Universität Erlangen-Nürnberg.
 *
 */
class IterationNumberBasedAdaptiveTimeStepping : public ITimeStepAlgorithm
{
public:

    /**
     * Constructor
     *
     * @param t_initial             start time
     * @param t_end                 end time
     * @param min_ts                the minimum allowed time step size
     * @param max_ts                the maximum allowed time step size
     * @param initial_ts            initial time step size
     * @param iter_times_vector     a vector of iteration numbers
     * (\f$i_1\f$, \f$i_2\f$, ..., \f$i_n\f$) which defines intervals as
     * \f$[i_1,i_2)\f$, \f$[i_2,i_3)\f$, ..., \f$[i_n,\infty)\f$.
     * If an iteration number is larger than \f$i_n\f$, current time step is
     * repeated with the new time step size.
     * @param multiplier_vector     a vector of multiplier coefficients
     * (\f$a_1\f$, \f$a_2\f$, ..., \f$a_n\f$) corresponding to the intervals given by iter_times_vector.
     * A time step size is calculated by \f$\Delta t_{n+1} = a * \Delta t_{n}\f$
     */
    IterationNumberBasedAdaptiveTimeStepping( double t_initial,
                                double t_end,
                                double min_ts,
                                double max_ts,
                                double initial_ts,
                                const std::vector<std::size_t> &iter_times_vector,
                                const std::vector<double> &multiplier_vector);

    virtual ~IterationNumberBasedAdaptiveTimeStepping() {}

    /// return the beginning of time steps
    virtual double begin() const {return _t_initial;}

    /// return the end of time steps
    virtual double end() const {return _t_end;}

    /// return current time step
    virtual const TimeStep getTimeStep() const;

    /// move to the next time step
    virtual bool next();

    /// return if the current step is accepted
    virtual bool accepted() const;

    /// return a history of time step sizes
    virtual const std::vector<double>& getTimeStepSizeHistory() const {return this->_dt_vector;}

    /// set the number of iterations
    void setNIterations(std::size_t n_itr) {this->_iter_times = n_itr;}

    /// return the number of repeated steps
    std::size_t getNumberOfRepeatedSteps() const {return this->_n_rejected_steps;}

private:
    /// calculate the next time step size
    double getNextTimeStepSize() const;

    /// initial time
    const double _t_initial;
    /// end time
    const double _t_end;
    /// this vector stores the number of iterations to which the respective multiplier coefficient will be applied
    const std::vector<std::size_t> _iter_times_vector;
    /// this vector stores the multiplier coefficients
    const std::vector<double> _multiplier_vector;
    /// the minimum allowed time step size
    const double _min_ts;
    /// the maximum allowed time step size
    const double _max_ts;
    /// initial time step size
    const double _initial_ts;
    /// the maximum allowed iteration number to accept current time step
    const std::size_t _max_iter;
    /// the number of nonlinear iterations
    std::size_t _iter_times;
    /// previous time step
    TimeStep _ts_pre;
    /// current time step
    TimeStep _ts_current;
    /// history of time step sizes
    std::vector<double> _dt_vector;
    /// the number of rejected steps
    std::size_t _n_rejected_steps;
};

} // NumLib
