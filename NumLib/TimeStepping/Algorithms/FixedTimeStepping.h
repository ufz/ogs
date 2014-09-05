/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef FIXEDTIMESTEPPING_H_
#define FIXEDTIMESTEPPING_H_

#include <vector>

#include "ITimeStepAlgorithm.h"

namespace NumLib
{

/**
 * \brief Fixed time stepping algorithm
 *
 * This algorithm returns time step size defined by a user priori.
 */
class FixedTimeStepping : public ITimeStepAlgorithm
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

    virtual ~FixedTimeStepping() {}

    /// return the beginning of time steps
    virtual double begin() const {return _t_initial;}

    /// return the end of time steps
    virtual double end() const {return _t_end;}

    /// return current time step
    virtual const TimeStep getTimeStep() const;

    /// move to the next time step
    virtual bool next();

    /// return if current time step is accepted
    virtual bool accepted() const {return true;}

    /// return a history of time step sizes
    virtual const std::vector<double>& getTimeStepSizeHistory() const {return _dt_vector; }

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
