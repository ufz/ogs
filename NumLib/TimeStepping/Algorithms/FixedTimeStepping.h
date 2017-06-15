/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

#include "TimeStepAlgorithm.h"

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
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
     *  n=\frac{t_{\rm end} - t_{\rm initial}}{\Delta t}
     * \f].
     *
     * @param t_initial     start time
     * @param t_end         finish time
     * @param dt            uniform time step size
     */
    FixedTimeStepping(double t_initial, double t_end, double dt);

    /**
     * Constructor with user-specified time step sizes
     *
     * A user can specify \f$\Delta t\f$ for each time step (i.e. \f$\Delta t_1,
     * \Delta t_2, ..., \Delta t_n\f$). Time at \f$m\f$ th step is given as
     * \f[
     *  t_{m}=\sum_{i=1}^m \Delta t_i + t_{\rm initial}
     * \f].
     *
     * @param t_initial     start time
     * @param t_end         finish time
     * @param vec_all_dt    a vector of all time steps
     */
    FixedTimeStepping(double t_initial, double t_end,
                      const std::vector<double>& vec_all_dt);

    /// Create timestepper from the given configuration
    static std::unique_ptr<TimeStepAlgorithm> newInstance(
        BaseLib::ConfigTree const& config);

    /// move to the next time step
    bool next(const double solution_error) override;

    /// return if current time step is accepted
    bool accepted() const override { return true; }
private:
    /// determine true end time
    static double computeEnd(double t_initial, double t_end,
                             const std::vector<double>& dt_vector);
};

}  // NumLib
