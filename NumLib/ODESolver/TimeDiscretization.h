/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "Types.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! Interface of time discretization schemes for first-order ODEs.
 *
 * The purpose of TimeDiscretization instances is to store the solution history
 * of an ODE, i. e., to keep the solution at as many timestamps as is required
 * by the respective time discretization scheme. Furthermore, TimeDiscretization
 * instances compute the discretized approximation of the time derivative \f$
 * \partial x/\partial t \f$.
 *
 * \note The method documentation of this class uses quantities introduced in
 * the following section.
 *
 * \todo Currently this interface does not yet support adaptive timestepping.
 *       While implementing that will lead to no changes for single-step
 * methods, for multi-step methods this interface will have to be extended.
 *
 *
 * Discretizing first-order ODEs {#concept_time_discretization}
 * =============================
 *
 * A first-order (implicit) ODE has the general form
 *
 * \f[ F(\dot x, x, t) \stackrel{!}{=} 0. \f]
 *
 * In order to solve it numerically a certain time discretization scheme,
 * such as the forward or backward Euler methods, is used.
 * The discretized ODE is then given by
 *
 * \f[ F(\hat x, x_C, t_C) \stackrel{!}{=} 0. \f]
 *
 * This interface has been designed with first-order implicit quasi-linear ODEs
 * in mind. They can be expressed as follows and are given here only to serve as
 * an example.
 *
 * \f[ M(x,t)\cdot \dot x + K(x,t) \cdot x - b(x,t)
 *  =: r(\dot x, x, t) \stackrel{!}{=} 0. \f]
 *
 * After time discretization this formula becomes:
 *
 * \f[ M(x_C,t_C)\cdot \hat x + K(x_C,t_C) \cdot x_C - b(x_C,t_C)
 *  =: r(\hat x, x_C, t_C) \stackrel{!}{=} 0. \f]
 *
 * The meaning of indices for \f$ x \f$ and \f$ t \f$ is as follows:
 *   * \f$ C \f$ -- "Current": The values of \f$ x \f$ and \f$ t \f$ at which
 * the discretized ODE is being assembled.
 *   * \f$ N \f$ -- "New": The values of \f$ x \f$ and \f$ t \f$ at the new
 * timestep that is being calculated right now by the ODE solver.
 *   * \f$ O \f$ -- "Old": The results from the preceding time step (or a linear
 * combination of results of the preceding time steps in the case of multistep
 * methods) weighted by a scalar factor.
 *   * \f$ n \f$ -- Numerical index indicating the timestep.
 *
 * \f$ \hat x \f$ is the discrete approximation of \f$ \dot x := \partial
 * x/\partial t\f$. It is assumed that \f$ \hat x \f$ can be written in the
 * following form: \f[ \hat x = \alpha \cdot x_N - x_O, \f] where \f$ \alpha :=
 * \partial \hat x / \partial x_N \f$ is a scalar.
 *
 * For different time discretization schemes \f$ x_C \f$, \f$ t_C \f$, \f$ x_N
 * \f$, \f$ x_O \f$ and \f$ \alpha \f$ take different values. Those for the time
 * implemented schemes are given in the table below.
 *
 * Scheme         | \f$ x_C \f$     | \f$ t_C \f$     | \f$ \alpha \f$        |
 * \f$ x_N \f$     | \f$ x_O \f$
 * -------------- | --------------- | --------------- | --------------------- |
 * --------------- | ---------------------- Backward Euler | \f$ x_{n+1} \f$ |
 * \f$ t_{n+1} \f$ | \f$ 1/\Delta t \f$    | \f$ x_{n+1} \f$ | \f$ x_n / \Delta
 * t \f$
 *
 */
class TimeDiscretization
{
public:
    TimeDiscretization() = default;

    //! Sets the initial condition.
    virtual void setInitialState(const double t0) = 0;

    /**
     * Compute and return the relative change of solutions between two
     * successive time steps by \f$ e_n = \|u^{n+1}-u^{n}\|/\|u^{n+1}\| \f$.
     *
     * @param x         The current solution
     * @param x_old     The previous solution
     * @param norm_type The norm type of global vector
     * @return          \f$ e_n = \|u^{n+1}-u^{n}\|/\|u^{n+1}\| \f$.
     *
     * @warning the value of x_old is changed to x - x_old after this
     * computation.
     */
    double computeRelativeChangeFromPreviousTimestep(
        GlobalVector const& x,
        GlobalVector const& x_old,
        MathLib::VecNormType norm_type);

    /*! Indicate that the computation of a new timestep is being started now.
     *
     * \warning Currently changing timestep sizes are not supported. Thus,
     *          \p delta_t must not change throughout the entire time
     *          integration process! This is not checked by this code!
     */
    virtual void nextTimestep(const double t, const double delta_t) = 0;

    //! Returns \f$ t_C \f$, i.e., the time at which the equation will be
    //! assembled.
    virtual double getCurrentTime() const = 0;

    //! Returns \f$ \Delta t_C \f$, i.e., the time at which the equation will be
    //! assembled.
    virtual double getCurrentTimeIncrement() const = 0;

    //! Returns \f$ \hat x \f$, i.e. the discretized approximation of \f$ \dot x
    //! \f$.
    void getXdot(GlobalVector const& x_at_new_timestep,
                 GlobalVector const& x_old,
                 GlobalVector& xdot) const
    {
        namespace LinAlg = MathLib::LinAlg;

        double const dt = getCurrentTimeIncrement();

        // xdot = 1/dt * x_at_new_timestep - x_old
        getWeightedOldX(xdot, x_old);
        LinAlg::axpby(xdot, 1. / dt, -1.0, x_at_new_timestep);
    }

    //! Returns \f$ x_O \f$.
    virtual void getWeightedOldX(
        GlobalVector& y, GlobalVector const& x_old) const = 0;  // = x_old

    virtual ~TimeDiscretization() = default;

protected:
    std::unique_ptr<GlobalVector> _dx;  ///< Used to store \f$ u_{n+1}-u_{n}\f$.
};

//! Backward Euler scheme.
class BackwardEuler final : public TimeDiscretization
{
public:
    void setInitialState(const double t0) override { _t = t0; }

    void nextTimestep(const double t, const double delta_t) override
    {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override { return _t; }
    double getCurrentTimeIncrement() const override { return _delta_t; }
    void getWeightedOldX(GlobalVector& y,
                         GlobalVector const& x_old) const override
    {
        namespace LinAlg = MathLib::LinAlg;

        // y = x_old / delta_t
        LinAlg::copy(x_old, y);
        LinAlg::scale(y, 1.0 / _delta_t);
    }

private:
    double _t = std::numeric_limits<double>::quiet_NaN();  //!< \f$ t_C \f$
    double _delta_t =
        std::numeric_limits<double>::quiet_NaN();  //!< the timestep size
};

//! @}
}  // namespace NumLib
