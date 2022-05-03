/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cmath>
#include <tuple>
#include <vector>

#include "BaseLib/Error.h"
#include "NumLib/TimeStepping/TimeStep.h"

namespace NumLib
{
/**
 * \brief Interface of time stepping algorithms
 */
class TimeStepAlgorithm
{
public:
    TimeStepAlgorithm(const double t0, const double t_end)
        : _t_initial(t0), _t_end(t_end), _ts_prev(t0), _ts_current(t0)
    {
    }

    virtual ~TimeStepAlgorithm() = default;

    /// return the beginning of time steps
    double begin() const { return _t_initial; }
    /// return the end of time steps
    double end() const { return _t_end; }
    /// return current time step
    const TimeStep getTimeStep() const { return _ts_current; }
    /// reset the current step size from the previous time
    virtual void resetCurrentTimeStep(const double dt)
    {
        updateTimesteps(dt, _ts_current, _ts_prev);
    }

    /// Move to the next time step
    /// \param solution_error Solution error \f$e_n\f$ between two successive
    ///        time steps.
    /// \param number_iterations Number of non-linear iterations used.
    /// \return A step acceptance flag and the computed step size.
    virtual std::tuple<bool, double> next(const double solution_error,
                                          int number_iterations) = 0;

    /// return if current time step is accepted or not
    virtual bool accepted() const { return _ts_current.isAccepted(); }

    /// Set the status of the step.
    /// \param accepted A boolean parameter is needed to indicated whether the
    /// step is accepted or not.
    void setAccepted(const bool accepted) { _ts_current.setAccepted(accepted); }

    /// Get a flag to indicate whether this algorithm needs to compute
    /// solution error. The default return value is false.
    virtual bool isSolutionErrorComputationNeeded() const { return false; }

    /// Query the timestepper if further time step size reduction is possible.
    virtual bool canReduceTimestepSize() const { return false; }

protected:
    /// initial time
    const double _t_initial;
    /// end time
    const double _t_end;

    /// previous time step information
    TimeStep _ts_prev;
    /// current time step information
    TimeStep _ts_current;
};

/// If any of the fixed times will be reached with given time increment, it will
/// be reduced, otherwise the input will be returned.
/// \pre The input vector of fixed times must be sorted.
/// \param t Current time.
/// \param dt Suggested time increment.
/// \param fixed_output_times Sorted list of times which are to be reached.
double possiblyClampDtToNextFixedTime(
    double const t, double const dt,
    std::vector<double> const& fixed_output_times);

}  // namespace NumLib
