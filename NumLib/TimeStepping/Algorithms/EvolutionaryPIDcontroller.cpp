/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on March 31, 2017, 4:13 PM
 */

#include "EvolutionaryPIDcontroller.h"

#include <functional>
#include <limits>
#include <vector>
#include "BaseLib/Logging.h"

#include "BaseLib/Algorithm.h"

namespace NumLib
{
EvolutionaryPIDcontroller::EvolutionaryPIDcontroller(
    const double t0, const double t_end, const double h0, const double h_min,
    const double h_max, const double rel_h_min, const double rel_h_max,
    std::vector<double>&& fixed_output_times, const double tol)
    : TimeStepAlgorithm(t0, t_end),
      h0_(h0),
      h_min_(h_min),
      h_max_(h_max),
      rel_h_min_(rel_h_min),
      rel_h_max_(rel_h_max),
      fixed_output_times_(std::move(fixed_output_times)),
      tol_(tol),
      e_n_minus1_(0.),
      e_n_minus2_(0.),
      is_accepted_(true)
{
    // Remove possible duplicated elements. Result will be sorted.
    BaseLib::makeVectorUnique(fixed_output_times_);
}

bool EvolutionaryPIDcontroller::next(double const solution_error,
                                     int const /*number_iterations*/)
{
    const bool is_previous_step_accepted = is_accepted_;

    const double e_n = solution_error;
    const double zero_threshlod = std::numeric_limits<double>::epsilon();
    // step rejected.
    if (e_n > tol_)
    {
        is_accepted_ = false;

        double h_new = (e_n > zero_threshlod) ? ts_current_.dt() * tol_ / e_n
                                              : 0.5 * ts_current_.dt();

        h_new = limitStepSize(h_new, is_previous_step_accepted);
        h_new = possiblyClampDtToNextFixedTime(ts_current_.current(), h_new,
                                               fixed_output_times_);

        ts_current_ = ts_prev_;
        ts_current_ += h_new;

        WARN(
            "This step is rejected due to the relative change from the"
            " solution of the previous\n"
            "\t time step to the current solution exceeds the given tolerance"
            " of {:g}.\n"
            "\t This time step will be repeated with a new time step size of"
            " {:g}\n"
            "\t or the simulation will be halted.",
            tol_, h_new);

        return false;
    }

    // step accepted.
    is_accepted_ = true;

    if (ts_current_.steps() == 0)
    {
        ts_prev_ = ts_current_;
        ts_current_ += h0_;
        e_n_minus1_ = e_n;

        dt_vector_.push_back(h0_);
    }
    else
    {
        const double h_n = ts_current_.dt();
        double h_new = h_n;

        if (e_n > zero_threshlod)
        {
            if (e_n_minus1_ > zero_threshlod)
            {
                if (e_n_minus2_ > zero_threshlod)
                {
                    h_new = std::pow(e_n_minus1_ / e_n, kP_) *
                            std::pow(tol_ / e_n, kI_) *
                            std::pow(
                                e_n_minus1_ * e_n_minus1_ / (e_n * e_n_minus2_),
                                kD_) *
                            h_n;
                }
                else
                {
                    h_new = std::pow(e_n_minus1_ / e_n, kP_) *
                            std::pow(tol_ / e_n, kI_) * h_n;
                }
            }
            else
            {
                h_new = std::pow(tol_ / e_n, kI_) * h_n;
            }
        }

        h_new = limitStepSize(h_new, is_previous_step_accepted);
        h_new = possiblyClampDtToNextFixedTime(ts_current_.current(), h_new,
                                               fixed_output_times_);
        dt_vector_.push_back(h_new);

        ts_prev_ = ts_current_;
        ts_current_ += h_new;

        e_n_minus2_ = e_n_minus1_;
        e_n_minus1_ = e_n;
    }

    return true;
}

double EvolutionaryPIDcontroller::limitStepSize(
    const double h_new, const bool previous_step_accepted) const
{
    const double h_n = ts_current_.dt();
    // Forced the computed time step size in the given range
    // (see the formulas in the documentation of the class)
    const double h_in_range = std::max(h_min_, std::min(h_new, h_max_));
    // Limit the step size change ratio.
    double limited_h =
        std::max(rel_h_min_ * h_n, std::min(h_in_range, rel_h_max_ * h_n));

    if (!previous_step_accepted)
    {
        // If the last time step was rejected and the new predicted time step
        // size is identical to that of the previous rejected step, the new
        // step size is then reduced by half.
        if (std::fabs(limited_h - ts_current_.dt()) <
            std::numeric_limits<double>::min())
        {
            limited_h = std::max(h_min_, 0.5 * limited_h);
        }

        // If the last time step was rejected and the new predicted time step
        // size is larger than the step size of the rejected step, the new step
        // size takes the half of the size of the rejected step. This could
        // happen when a time step is rejected due to a diverged non-linear
        // solver. In such case, this algorithm may give a large time step size
        // by using the diverged solution.
        if (limited_h > ts_current_.dt())
        {
            limited_h = std::max(h_min_, 0.5 * ts_current_.dt());
        }
    }
    return limited_h;
}

void EvolutionaryPIDcontroller::addFixedOutputTimes(
    std::vector<double> const& extra_fixed_output_times)
{
    fixed_output_times_.insert(fixed_output_times_.end(),
                               extra_fixed_output_times.begin(),
                               extra_fixed_output_times.end());

    // Remove possible duplicated elements. Result will be sorted.
    BaseLib::makeVectorUnique(fixed_output_times_);
}

bool EvolutionaryPIDcontroller::canReduceTimestepSize() const
{
    // If current and previous dt are both at minimum dt, then cannot reduce
    // further.
    return !(ts_current_.dt() == h_min_ && ts_prev_.dt() == h_min_);
}
}  // namespace NumLib
