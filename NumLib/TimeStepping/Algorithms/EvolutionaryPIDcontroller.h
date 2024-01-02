/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 31, 2017, 4:13 PM
 */

#pragma once

#include <memory>
#include <vector>

#include "TimeStepAlgorithm.h"

namespace NumLib
{
/**
 *  This class gives an adaptive algorithm whose time step control is
 *  evolutionary PID controller. With an definition of relative solution change
 *  \f$e_n=\frac{\|u^{n+1}-u^n\|}{\|u^{n+1}\|}\f$, the algorithm gives a time
 *   step size estimation as
 *   \f[
 *     h_{n+1} =  \left(\frac{e_{n-1}}{e_n}\right)^{k_P}
 *                \left(\frac{TOL}{e_n}\right)^{k_I}
 *                \left(\frac{e^2_{n-1}}{e_n e_{n-2}}\right)^{k_D}
 *   \f]
 *   where \f$k_P=0.075\f$, \f$k_I=0.175\f$, \f$k_D=0.01\f$ are empirical PID
 *   parameters.
 *
 *   In the computation, \f$ e_n\f$ is calculated firstly. If \f$e_n>TOL\f$, the
 *   current time step is rejected and repeated with a new time step size of
 *   \f$h=\frac{TOL}{e_n} h_n\f$.
 *
 *   Limits of the time step size are given as
 *   \f[
 *        h_{\mbox{min}} \leq h_{n+1} \leq h_{\mbox{max}},
 *        l \leq \frac{h_{n+1}}{h_n} \leq L
 *   \f]
 *
 *   Similar algorithm can be found in \cite ahmed2015adaptive .
 */
class EvolutionaryPIDcontroller final : public TimeStepAlgorithm
{
public:
    EvolutionaryPIDcontroller(const double t0, const double t_end,
                              const double h0, const double h_min,
                              const double h_max, const double rel_h_min,
                              const double rel_h_max, const double tol,
                              std::vector<double> const& fixed_times_for_output)
        : TimeStepAlgorithm(t0, t_end),
          _h0(h0),
          _h_min(h_min),
          _h_max(h_max),
          _rel_h_min(rel_h_min),
          _rel_h_max(rel_h_max),
          _tol(tol),
          _e_n_minus1(0.),
          _e_n_minus2(0.),
          _fixed_times_for_output(fixed_times_for_output)
    {
    }

    std::tuple<bool, double> next(double solution_error,
                                  int number_iterations,
                                  NumLib::TimeStep& timestep_previous,
                                  NumLib::TimeStep& timestep_current) override;

    bool isSolutionErrorComputationNeeded() const override { return true; }

    virtual bool canReduceTimestepSize(
        NumLib::TimeStep const& timestep_previous,
        NumLib::TimeStep const& timestep_current) const override;

private:
    const double _kP = 0.075;  ///< Parameter. \see EvolutionaryPIDcontroller
    const double _kI = 0.175;  ///< Parameter. \see EvolutionaryPIDcontroller
    const double _kD = 0.01;   ///< Parameter. \see EvolutionaryPIDcontroller

    const double _h0;     ///< initial time step size.
    const double _h_min;  ///< minimum step size.
    const double _h_max;  ///< maximum step size.

    /// \f$l\f$ in \f$ h_{\mbox{min}} \leq h_{n+1} \leq h_{\mbox{max}},\f$
    const double _rel_h_min;
    /// \f$L\f$ in \f$ h_{\mbox{min}} \leq h_{n+1} \leq h_{\mbox{max}},\f$
    const double _rel_h_max;

    const double _tol;

    double _e_n_minus1;  ///< \f$e_{n-1}\f$.
    double _e_n_minus2;  ///< \f$e_{n-2}\f$.

    std::vector<double> const _fixed_times_for_output;

    /**
     * Force the computed time step size in the given range
     * (see the formulas in the documentation of the class)
     * or use the half of the previous time step size under some other
     * constrains.
     * @param h_new                   The computed time step size.
     * @param previous_step_accepted  An indicator for whether the previous time
     *                                step is rejected.
     * @param timestep_current the current time step
     * @return                        The new time step after apply
     *                                the constrains.
     */
    double limitStepSize(const double h_new,
                         const bool previous_step_accepted,
                         NumLib::TimeStep const& timestep_current) const;
};

}  // end of namespace NumLib
