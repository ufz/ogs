/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
        : _t_initial(t0), _t_end(t_end)
    {
    }

    virtual ~TimeStepAlgorithm() = default;

    /// return the beginning of time steps
    double begin() const { return _t_initial; }
    /// return the end of time steps
    double end() const { return _t_end; }
    /// reset the current step size from the previous time
    virtual void resetCurrentTimeStep(const double /*dt*/,
                                      TimeStep& /*ts_previous*/,
                                      TimeStep& /*ts_current*/)
    {
    }

    /// Move to the next time step
    /// \param solution_error Solution error \f$e_n\f$ between two successive
    ///        time steps.
    /// \param number_iterations Number of non-linear iterations used.
    /// \param ts_previous the previous time step used to compute the size of
    /// the next step
    /// \param ts_current the current time step used to compute the size of the
    /// next step
    /// \return A step acceptance flag and the computed step size.
    virtual std::tuple<bool, double> next(const double solution_error,
                                          int number_iterations,
                                          NumLib::TimeStep& ts_previous,
                                          NumLib::TimeStep& ts_current) = 0;

    /// Get a flag to indicate whether this algorithm needs to compute
    /// solution error. The default return value is false.
    virtual bool isSolutionErrorComputationNeeded() const { return false; }

    /// Query the timestepper if further time step size reduction is possible.
    virtual bool canReduceTimestepSize(
        NumLib::TimeStep const& /*timestep_previous*/,
        NumLib::TimeStep const& /*timestep_current*/) const
    {
        return false;
    }

protected:
    /// initial time
    const double _t_initial;
    /// end time
    const double _t_end;
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

bool canReduceTimestepSize(TimeStep const& timestep_previous,
                           TimeStep const& timestep_current,
                           double const min_dt);
}  // namespace NumLib
