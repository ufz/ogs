/**
 * \file
 * \author Haibing Shao and Norihiro Watanabe
 * \date   2013-08-07
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IterationNumberBasedTimeStepping.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

#include "BaseLib/Algorithm.h"

namespace NumLib
{
IterationNumberBasedTimeStepping::IterationNumberBasedTimeStepping(
    double const t_initial, double const t_end, double const min_dt,
    double const max_dt, double const initial_dt,
    std::vector<int>&& iter_times_vector,
    std::vector<double>&& multiplier_vector,
    std::vector<double>&& fixed_output_times)
    : TimeStepAlgorithm(t_initial, t_end),
      iter_times_vector_(std::move(iter_times_vector)),
      multiplier_vector_(std::move(multiplier_vector)),
      fixed_output_times_(std::move(fixed_output_times)),
      min_dt_(min_dt),
      max_dt_(max_dt),
      initial_dt_(initial_dt),
      max_iter_(iter_times_vector_.empty() ? 0 : iter_times_vector_.back())
{
    if (iter_times_vector_.empty())
    {
        OGS_FATAL("Vector of iteration numbers must not be empty.");
    }
    if (iter_times_vector_.size() != multiplier_vector_.size())
    {
        OGS_FATAL(
            "Vector of iteration numbers must be of the same size as the "
            "vector of multipliers.");
    }
    if (!std::is_sorted(std::begin(iter_times_vector_),
                        std::end(iter_times_vector_)))
    {
        OGS_FATAL("Vector of iteration numbers must be sorted.");
    }

    // Remove possible duplicated elements. Result will be sorted.
    BaseLib::makeVectorUnique(fixed_output_times_);
}

bool IterationNumberBasedTimeStepping::next(double const /*solution_error*/,
                                            int const number_iterations)
{
    iter_times_ = number_iterations;

    // confirm current time and move to the next if accepted
    if (accepted())
    {
        ts_prev_ = ts_current_;
        dt_vector_.push_back(ts_current_.dt());
    }
    else
    {
        ++n_rejected_steps_;
        // time step was rejected, keep dt for the next dt computation.
        ts_prev_ =  // essentially equal to ts_prev_.dt = ts_current_.dt.
            TimeStep{ts_prev_.previous(),
                     ts_prev_.previous() + ts_current_.dt(), ts_prev_.steps()};
    }

    // prepare the next time step info
    ts_current_ = ts_prev_;
    ts_current_ += possiblyClampDtToNextFixedTime(
        ts_current_.current(), getNextTimeStepSize(), fixed_output_times_);

    return true;
}

double IterationNumberBasedTimeStepping::findMultiplier(
    int const number_iterations) const
{
    double multiplier = multiplier_vector_.front();
    for (std::size_t i = 0; i < iter_times_vector_.size(); i++)
    {
        if (number_iterations >= iter_times_vector_[i])
        {
            multiplier = multiplier_vector_[i];
        }
    }

    if (!accepted_ && (multiplier >= 1.0))
    {
        return *std::min_element(multiplier_vector_.begin(),
                                 multiplier_vector_.end());
    }

    return multiplier;
}

double IterationNumberBasedTimeStepping::getNextTimeStepSize() const
{
    double dt = 0.0;

    // In first time step and first non-linear iteration take the initial dt.
    if (ts_prev_.steps() == 0 && iter_times_ == 0)
    {
        dt = initial_dt_;
    }
    else
    {
        // Attention: for the first time step and second iteration the
        // ts_prev.dt is 0 and 0*multiplier is the next dt, which will be
        // clamped to the minimum dt.
        dt = ts_prev_.dt() * findMultiplier(iter_times_);
    }

    dt = std::clamp(dt, min_dt_, max_dt_);

    double const t_next = dt + ts_prev_.current();
    if (t_next > end())
    {
        dt = end() - ts_prev_.current();
    }

    return dt;
}

void IterationNumberBasedTimeStepping::addFixedOutputTimes(
    std::vector<double> const& extra_fixed_output_times)
{
    fixed_output_times_.insert(fixed_output_times_.end(),
                               extra_fixed_output_times.begin(),
                               extra_fixed_output_times.end());

    // Remove possible duplicated elements. Result will be sorted.
    BaseLib::makeVectorUnique(fixed_output_times_);
}

bool IterationNumberBasedTimeStepping::accepted() const
{
    return accepted_;
}

bool IterationNumberBasedTimeStepping::canReduceTimestepSize() const
{
    // If current and previous dt are both at minimum dt, then cannot reduce
    // further.
    return !(ts_current_.dt() == min_dt_ && ts_prev_.dt() == min_dt_);
}

}  // namespace NumLib
