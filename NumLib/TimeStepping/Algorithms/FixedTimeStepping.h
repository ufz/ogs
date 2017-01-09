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

#ifndef FIXEDTIMESTEPPING_H_
#define FIXEDTIMESTEPPING_H_

#include <memory>
#include <vector>

#include "ITimeStepAlgorithm.h"

namespace BaseLib { class ConfigTree; }

namespace NumLib
{

/**
 * \brief Fixed time stepping algorithm
 *
 * This algorithm returns time step size defined by a user priori.
 */
class FixedTimeStepping final
        : public ITimeStepAlgorithm
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
    FixedTimeStepping(double t_initial, double t_end, const std::vector<double> &vec_all_dt);

    /// Create timestepper from the given configuration
    static std::unique_ptr<ITimeStepAlgorithm> newInstance(BaseLib::ConfigTree const& config);

    /// return the beginning of time steps
    double begin() const override { return _t_initial; }

    /// return the end of time steps
    double end() const override { return _t_end; }

    /// return current time step
    const TimeStep getTimeStep() const override;

    /// move to the next time step
    bool next() override;

    /// return if current time step is accepted
    bool accepted() const override { return true; }

    /// return a history of time step sizes
    const std::vector<double>& getTimeStepSizeHistory() const override { return _dt_vector; }

private:
    /// determine true end time
    static double computeEnd(double t_initial, double t_end, const std::vector<double> &dt_vector);

    /// initial time
    const double _t_initial;
    /// end time
    const double _t_end;
    /// a vector of time step sizes
    const std::vector<double> _dt_vector;
    /// previous time step information
    TimeStep _ts_prev;
    /// current time step information
    TimeStep _ts_current;
};

} //NumLib

#endif // FIXEDTIMESTEPPING_H_
