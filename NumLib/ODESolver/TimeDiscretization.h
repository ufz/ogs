/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

//! Interface that allows managing an additional internal state as required by
//! certain time discretization schemes.
class InternalMatrixStorage
{
public:
    /*! Triggers a refresh of the internal matrix/vector storage.
     *
     * \remark
     * This method is needed in particular to fully implement the
     * interaction of the CrankNicolson scheme with other classes.
     *
     * \attention
     * This method must be called (if it is called) from within
     * TimeDiscretization::pushState() \b after the internal state of
     * the TimeDiscretization has been set to the new solution.
     * Otherwise the pushMatrices() method of MatrixTranslator's will break!
     */
    virtual void pushMatrices() const = 0;

    virtual ~InternalMatrixStorage() = default;
};

/*! Interface of time discretization schemes for first-order ODEs.
 *
 * The purpose of TimeDiscretization instances is to store the solution history of
 * an ODE, i. e., to keep the solution at as many timestamps as is required by the
 * respective time discretization scheme. Furthermore, TimeDiscretization instances
 * compute the discretized approximation of the time derivative \f$ \partial x/\partial t \f$.
 *
 * \note The method documentation of this class uses quantities introduced in the
 *       following section.
 *
 * \todo Currently this interface does not yet support adaptive timestepping.
 *       While implementing that will lead to no changes for single-step methods,
 *       for multi-step methods this interface will have to be extended.
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
 * This interface has been designed with first-order implicit quasi-linear ODEs in mind.
 * They can be expressed as follows and are given here only to serve as an example.
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
 *   * \f$ C \f$ -- "Current": The values of \f$ x \f$ and \f$ t \f$ at which the discretized
 *                             ODE is being assembled.
 *   * \f$ N \f$ -- "New": The values of \f$ x \f$ and \f$ t \f$ at the new timestep that is
 *                         being calculated right now by the ODE solver.
 *   * \f$ O \f$ -- "Old": The results from the preceding time step (or a linear combination of
 *                         results of the preceding time steps in the case of multistep methods)
 *                         weighted by a scalar factor.
 *   * \f$ n \f$ -- Numerical index indicating the timestep.
 *
 * \f$ \hat x \f$ is the discrete approximation of \f$ \dot x := \partial x/\partial t\f$.
 * It is assumed that \f$ \hat x \f$ can be written in the following form:
 * \f[ \hat x = \alpha \cdot x_N - x_O, \f]
 * where \f$ \alpha := \partial \hat x / \partial x_N \f$ is a scalar.
 *
 * For different time discretization schemes \f$ x_C \f$, \f$ t_C \f$, \f$ x_N \f$,
 * \f$ x_O \f$ and \f$ \alpha \f$ take different values.
 * Those for the time implemented schemes are given in the table below.
 *
 * Scheme         | \f$ x_C \f$     | \f$ t_C \f$     | \f$ \alpha \f$        | \f$ x_N \f$     | \f$ x_O \f$
 * -------------- | --------------- | --------------- | --------------------- | --------------- | ----------------------
 * Forward Euler  | \f$ x_n \f$     | \f$ t_n \f$     | \f$ 1/\Delta t \f$    | \f$ x_{n+1} \f$ | \f$ x_n / \Delta t \f$
 * Backward Euler | \f$ x_{n+1} \f$ | \f$ t_{n+1} \f$ | \f$ 1/\Delta t \f$    | \f$ x_{n+1} \f$ | \f$ x_n / \Delta t \f$
 * Crank-Nicolson | \f$ x_{n+1} \f$ | \f$ t_{n+1} \f$ | \f$ 1/\Delta t \f$    | \f$ x_{n+1} \f$ | \f$ x_n / \Delta t \f$
 * BDF(2)         | \f$ x_{n+2} \f$ | \f$ t_{n+1} \f$ | \f$ 3/(2\Delta t) \f$ | \f$ x_{n+2} \f$ | \f$ (2\cdot x_{n+1} - x_n/2)/\Delta t \f$
 *
 * The other backward differentiation formulas of orders 1 to 6 are also implemented, but only
 * BDF(2) has bee given here for brevity.
 *
 */
class TimeDiscretization
{
public:
    //! Sets the initial condition.
    virtual void setInitialState(const double t0, GlobalVector const& x0) = 0;

    /*! Indicate that the current timestep is done and that you will proceed to
     * the next one.
     *
     * \warning Do not use this method for setting the initial condition,
     *          rather use setInitialState()!
     *
     * \param t    The current timestep.
     * \param x    The solution at the current timestep.
     * \param strg Trigger storing some internal state.
     *             Currently only used by the CrankNicolson scheme.
     */
    virtual void pushState(const double t, GlobalVector const& x,
                           InternalMatrixStorage const& strg) = 0;

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

    //! Returns \f$ \hat x \f$, i.e. the discretized approximation of \f$ \dot x
    //! \f$.
    void getXdot(GlobalVector const& x_at_new_timestep, GlobalVector& xdot) const
    {
        namespace LinAlg = MathLib::LinAlg;

        auto const dxdot_dx = getNewXWeight();

        // xdot = dxdot_dx * x_at_new_timestep - x_old
        getWeightedOldX(xdot);
        LinAlg::axpby(xdot, dxdot_dx, -1.0, x_at_new_timestep);
    }

    //! Returns \f$ \alpha = \partial \hat x / \partial x_N \f$.
    virtual double getNewXWeight() const = 0;

    //! Returns \f$ x_O \f$.
    virtual void getWeightedOldX(GlobalVector& y) const = 0;  // = x_old

    virtual ~TimeDiscretization() = default;

    //! \name Extended Interface
    //! These methods are provided primarily to make certain concrete time
    //! discretizations
    //! with special demands, such as the forward Euler or Crank-Nicolson
    //! schemes, possible.
    //! @{

    /*! Tell whether this scheme inherently requires a nonlinear solver or not.
     *
     * The ForwardEuler scheme is inherently linear in that sense, the others
     * are not.
     */
    virtual bool isLinearTimeDisc() const { return false; }
    /*! Returns \f$ \partial x_C / \partial x_N \f$.
     *
     * The ForwardEuler scheme overrides this.
     */
    virtual double getDxDx() const { return 1.0; }
    /*! Returns \f$ x_C \f$, i.e., the state at which the equation will be
     * assembled.
     *
     * This method is overridden in the ForwardEuler scheme.
     */
    virtual GlobalVector const& getCurrentX(GlobalVector const& x_at_new_timestep) const
    {
        return x_at_new_timestep;
    }

    /*! Indicate that this scheme needs some additional assembly before the
     * first
     *  timestep will be solved.
     *
     * The CrankNicolson scheme needs such preload.
     */
    virtual bool needsPreload() const { return false; }
    //! @}
};

//! Backward Euler scheme.
class BackwardEuler final : public TimeDiscretization
{
public:
    BackwardEuler()
        : _x_old(NumLib::GlobalVectorProvider::provider.getVector())
    {
    }

    ~BackwardEuler()
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(_x_old);
    }

    void setInitialState(const double t0, GlobalVector const& x0) override
    {
        _t = t0;
        MathLib::LinAlg::copy(x0, _x_old);
    }

    void pushState(const double /*t*/, GlobalVector const& x,
                   InternalMatrixStorage const&) override
    {
        MathLib::LinAlg::copy(x, _x_old);
    }

    void nextTimestep(const double t, const double delta_t) override
    {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override { return _t; }
    double getNewXWeight() const override { return 1.0 / _delta_t; }
    void getWeightedOldX(GlobalVector& y) const override
    {
        namespace LinAlg = MathLib::LinAlg;

        // y = x_old / delta_t
        LinAlg::copy(_x_old, y);
        LinAlg::scale(y, 1.0 / _delta_t);
    }

private:
    double _t;        //!< \f$ t_C \f$
    double _delta_t;  //!< the timestep size
    GlobalVector& _x_old;   //!< the solution from the preceding timestep
};

//! Forward Euler scheme.
class ForwardEuler final : public TimeDiscretization
{
public:
    ForwardEuler()
        : _x_old(NumLib::GlobalVectorProvider::provider.getVector())
    {
    }

    ~ForwardEuler()
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(_x_old);
    }

    void setInitialState(const double t0, GlobalVector const& x0) override
    {
        _t = t0;
        _t_old = t0;
        MathLib::LinAlg::copy(x0, _x_old);
    }

    void pushState(const double /*t*/, GlobalVector const& x,
                   InternalMatrixStorage const&) override
    {
        MathLib::LinAlg::copy(x, _x_old);
    }

    void nextTimestep(const double t, const double delta_t) override
    {
        _t_old = _t;
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override
    {
        return _t_old;  // forward Euler does assembly at the preceding timestep
    }

    GlobalVector const& getCurrentX(
        const GlobalVector& /*x_at_new_timestep*/) const override
    {
        return _x_old;
    }

    double getNewXWeight() const override { return 1.0 / _delta_t; }
    void getWeightedOldX(GlobalVector& y) const override
    {
        namespace LinAlg = MathLib::LinAlg;

        // y = x_old / delta_t
        LinAlg::copy(_x_old, y);
        LinAlg::scale(y, 1.0 / _delta_t);
    }

    bool isLinearTimeDisc() const override { return true; }
    double getDxDx() const override { return 0.0; }
    //! Returns the solution from the preceding timestep.
    GlobalVector const& getXOld() const { return _x_old; }
private:
    double _t;        //!< \f$ t_C \f$
    double _t_old;    //!< the time of the preceding timestep
    double _delta_t;  //!< the timestep size
    GlobalVector& _x_old;   //!< the solution from the preceding timestep
};

//! Generalized Crank-Nicolson scheme.
class CrankNicolson final : public TimeDiscretization
{
public:
    /*! Constructs a new instance.
     *
     * \param theta The implicitness parameter \f$ \theta \f$. Some special
     * values are:
     *              \arg 1.0 fully implicit (like BackwardEuler).
     *              \arg 0.0 fully explicit (like ForwardEuler).
     *              \arg 0.5 traditional Crank-Nicolson scheme.
     */
    explicit CrankNicolson(const double theta)
        : _theta(theta),
          _x_old(NumLib::GlobalVectorProvider::provider.getVector())
    {
    }

    ~CrankNicolson()
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(_x_old);
    }

    void setInitialState(const double t0, GlobalVector const& x0) override
    {
        _t = t0;
        MathLib::LinAlg::copy(x0, _x_old);
    }

    void pushState(const double, GlobalVector const& x,
                   InternalMatrixStorage const& strg) override
    {
        MathLib::LinAlg::copy(x, _x_old);
        strg.pushMatrices();
    }

    void nextTimestep(const double t, const double delta_t) override
    {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override { return _t; }
    double getNewXWeight() const override { return 1.0 / _delta_t; }
    void getWeightedOldX(GlobalVector& y) const override
    {
        namespace LinAlg = MathLib::LinAlg;

        // y = x_old / delta_t
        LinAlg::copy(_x_old, y);
        LinAlg::scale(y, 1.0 / _delta_t);
    }

    bool needsPreload() const override { return true; }
    //! Returns \f$ \theta \f$.
    double getTheta() const { return _theta; }
    //! Returns the solution from the preceding timestep.
    GlobalVector const& getXOld() const { return _x_old; }
private:
    const double _theta;  //!< the implicitness parameter \f$ \theta \f$
    double _t;            //!< \f$ t_C \f$
    double _delta_t;      //!< the timestep size
    GlobalVector& _x_old;       //!< the solution from the preceding timestep
};

namespace detail
{
//! Coefficients used in the backward differentiation formulas.
const double BDF_Coeffs[6][7] = {
    // leftmost column: weight of the solution at the new timestep
    // signs of columns > 1 are flipped compared to standard BDF tableaus
    {1.0, 1.0},
    {1.5, 2.0, -0.5},
    {11.0 / 6.0, 3.0, -1.5, 1.0 / 3.0},
    {25.0 / 12.0, 4.0, -3.0, 4.0 / 3.0, -0.25},
    {137.0 / 60.0, 5.0, -5.0, 10.0 / 3.0, -1.25, 1.0 / 5.0},
    {147.0 / 60.0, 6.0, -7.5, 20.0 / 3.0, -3.75, 6.0 / 5.0, -1.0 / 6.0}
    // coefficient of (for BDF(6), the oldest state, x_n, is always rightmost)
    //        x_+6, x_+5, x_+4,       x_+3,  x_+2, x_+1,     x_n
};
}

//! Backward differentiation formula.
class BackwardDifferentiationFormula final : public TimeDiscretization
{
public:
    /*! Constructs a new instance.
     *
     * \param num_steps The order of the BDF to be used
     *                  (= the number of timesteps kept in the internal history
     * buffer).
     *                  Valid range: 1 through 6.
     *
     * \note Until a sufficient number of timesteps has been computed to be able
     *       to use the full \c num_steps order BDF, lower order BDFs are used
     * in
     *       the first timesteps.
     */
    explicit BackwardDifferentiationFormula(const unsigned num_steps)
        : _num_steps(num_steps)
    {
        assert(1 <= num_steps && num_steps <= 6);
        _xs_old.reserve(num_steps);
    }

    ~BackwardDifferentiationFormula()
    {
        for (auto* x : _xs_old)
            NumLib::GlobalVectorProvider::provider.releaseVector(*x);
    }

    void setInitialState(const double t0, GlobalVector const& x0) override
    {
        _t = t0;
        _xs_old.push_back(
            &NumLib::GlobalVectorProvider::provider.getVector(x0));
    }

    void pushState(const double, GlobalVector const& x,
                   InternalMatrixStorage const&) override
    {
        namespace LinAlg = MathLib::LinAlg;
        // TODO use boost cirular buffer?

        // until _xs_old is filled, lower-order BDF formulas are used.
        if (_xs_old.size() < _num_steps)
        {
            _xs_old.push_back(
                &NumLib::GlobalVectorProvider::provider.getVector(x));
        }
        else
        {
            LinAlg::copy(x, *_xs_old[_offset]);
            _offset = (_offset + 1) %
                      _num_steps;  // treat _xs_old as a circular buffer
        }
    }

    void nextTimestep(const double t, const double delta_t) override
    {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override { return _t; }
    double getNewXWeight() const override
    {
        auto const k = eff_num_steps();
        return detail::BDF_Coeffs[k - 1][0] / _delta_t;
    }

    void getWeightedOldX(GlobalVector& y) const override
    {
        namespace LinAlg = MathLib::LinAlg;

        auto const k = eff_num_steps();
        auto const* const BDFk = detail::BDF_Coeffs[k - 1];

        // compute linear combination \sum_{i=0}^{k-1} BDFk_{k-i} \cdot x_{n+i}
        LinAlg::copy(*_xs_old[_offset], y);  // _xs_old[offset] = x_n
        LinAlg::scale(y, BDFk[k]);

        for (unsigned i = 1; i < k; ++i)
        {
            auto const off = (_offset + i) % k;
            LinAlg::axpy(y, BDFk[k - i], *_xs_old[off]);
        }

        LinAlg::scale(y, 1.0 / _delta_t);
    }

private:
    std::size_t eff_num_steps() const { return _xs_old.size(); }
    const unsigned _num_steps;  //!< The order of the BDF method
    double _t;                  //!< \f$ t_C \f$
    double _delta_t;            //!< the timestep size

    std::vector<GlobalVector*> _xs_old;  //!< solutions from the preceding timesteps
    unsigned _offset = 0;  //!< allows treating \c _xs_old as circular buffer
};

//! @}
}
