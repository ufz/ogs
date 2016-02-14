#pragma once

#include <vector>

#include "MathLib/LinAlg/BLAS.h"
#include "Types.h"


namespace NumLib
{

//! \addtogroup ODESolver
//! @{

class InternalMatrixStorage
{
public:
    // needed for Crank-Nicolson
    virtual void pushMatrices() const = 0;
};

/*! Interface of time discretization schemes.
 *
 * This interface can be used to describe time discretization schemes for
 * first-order ODEs.
 *
 * The purpose of TimeDiscretization instances is to store the solution history of
 * an ODE, i. e., to keep the solution at as many timestamps as is required by the
 * respective time discretization scheme. Furthermore, TimeDiscretization instances
 * compute the discretized approximation of the time derivative \f$ \partial x/\partial t \f$.
 *
 *
 * Design Ideas
 * ------------
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
 * \f[ \hat x = \alpha \cdot x_C - x_O, \f]
 * where \f$ \alpha := \partial \hat x / \partial x_C \f$ is a scalar.
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
template<typename Vector>
class TimeDiscretization
{
public:
    virtual void setInitialState(const double t0, Vector const& x) = 0;
    virtual void pushState(const double t, Vector const& x,
                           InternalMatrixStorage const& eq) = 0;

    virtual void setCurrentTime(const double t, const double delta_t) = 0;
    virtual double getCurrentTime() const = 0; // get time used for assembly

    // \dot x === alpha * x - x_old
    void getXdot(Vector const& x_at_new_timestep, Vector& xdot) const
    {
        namespace BLAS = MathLib::BLAS;

        auto const dxdot_dx = getCurrentXWeight();

        getWeightedOldX(xdot);
        BLAS::axpby(xdot, dxdot_dx, -1.0, x_at_new_timestep);
    }

    virtual double getCurrentXWeight() const = 0; // = alpha
    virtual void getWeightedOldX(Vector& y) const = 0; // = x_old

    ~TimeDiscretization() = default;

    //! \name Extended Interface
    //! These methods are provided primarily to make certain concrete time discretizations
    //! with special demands, such as the forward Euler or Crank-Nicolson schemes, possible.
    //! @{

    // Forward Euler is linear, other schemes not.
    virtual bool isLinearTimeDisc() const { return false; }

    // Forward Euler will override this
    virtual double getDxDx() const { return 1.0; }

    // Forward Euler overrides this.
    // Caution: This is not the x with which you want to compute \dot x
    virtual Vector const& getCurrentX(Vector const& x_at_new_timestep) const
    {
        return x_at_new_timestep;
    }

    // for Crank-Nicolson
    virtual bool needsPreload() const { return false; }

    //! @}
};


template<typename Vector>
class BackwardEuler final : public TimeDiscretization<Vector>
{
public:
    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _x_old = x;
    }

    void pushState(const double t, Vector const& x, InternalMatrixStorage const&) override
    {
        (void) t;
        _x_old = x;
    }

    void setCurrentTime(const double t, const double delta_t) override {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override {
        return _t;
    }

    double getCurrentXWeight() const override {
        return 1.0/_delta_t;
    }

    void getWeightedOldX(Vector& y) const override
    {
        namespace BLAS = MathLib::BLAS;
        BLAS::copy(_x_old, y);
        BLAS::scale(y, 1.0/_delta_t);
    }

private:
    double _t = 9999.9999;
    double _delta_t = 8888.8888;
    Vector _x_old;
};


template<typename Vector>
class ForwardEuler final : public TimeDiscretization<Vector>
{
public:
    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _t_old = t0;
        _x_old = x;
    }

    void pushState(const double t, Vector const& x, InternalMatrixStorage const&) override
    {
        (void) t;
        _x_old = x;
    }

    void setCurrentTime(const double t, const double delta_t) override {
        _t_old = _t;
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override {
        return _t_old; // forward Euler does assembly at the preceding timestep
    }

    Vector const& getCurrentX(const Vector& /*x_at_new_timestep*/) const override {
        return _x_old;
    }

    double getCurrentXWeight() const override {
        return 1.0/_delta_t;
    }

    void getWeightedOldX(Vector& y) const override
    {
        namespace BLAS = MathLib::BLAS;
        BLAS::copy(_x_old, y);
        BLAS::scale(y, 1.0/_delta_t);
    }

    bool isLinearTimeDisc() const override {
        return true;
    }

    double getDxDx() const override {
        return 0.0;
    }

    Vector const& getXOld() const { return _x_old; }

private:
    double _t = 9999.9999;
    double _t_old = 7777.7777;
    double _delta_t = 8888.8888;
    Vector _x_old;
};


template<typename Vector>
class CrankNicolson final : public TimeDiscretization<Vector>
{
public:
    explicit
    CrankNicolson(const double theta)
        : _theta(theta)
    {}

    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _x_old = x;
    }

    void pushState(const double t, Vector const& x, InternalMatrixStorage const& eq) override
    {
        (void) t;
        _x_old = x;
        eq.pushMatrices();
    }

    void setCurrentTime(const double t, const double delta_t) override {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override {
        return _t;
    }

    double getCurrentXWeight() const override {
        return 1.0/_delta_t;
    }

    void getWeightedOldX(Vector& y) const override
    {
        namespace BLAS = MathLib::BLAS;
        BLAS::copy(_x_old, y);
        BLAS::scale(y, 1.0/_delta_t);
    }

    bool needsPreload() const override {
        return true;
    }

    double getTheta() const { return _theta; }
    Vector const& getXOld() const { return _x_old; }

private:
    const double _theta = 555.555;
    double _t = 9999.9999;
    double _delta_t = 8888.8888;
    Vector _x_old;
};


namespace detail
{

const double BDF_Coeffs[6][7] = {
    // leftmost column: weight of the solution at the new timestep
    // signs of columns > 1 are flipped compared to standard BDF tableaus
    {   1.0,         1.0 },
    {   1.5,         2.0, -0.5 },
    {  11.0 /  6.0,  3.0, -1.5,  1.0 / 3.0 },
    {  25.0 / 12.0,  4.0, -3.0,  4.0 / 3.0, -0.25 },
    { 137.0 / 60.0,  5.0, -5.0, 10.0 / 3.0, -1.25,  0.2 },
    { 147.0 / 60.0,  6.0, -7.5, 20.0 / 3.0, -3.75,  1.2, -1.0/6.0 }
    // coefficient of (for BDF(6), the oldest state, x_n, is always rightmost)
    //        x_+6, x_+5, x_+4,       x_+3,  x_+2, x_+1,     x_n
};

}


template<typename Vector>
class BackwardDifferentiationFormula final : public TimeDiscretization<Vector>
{
public:
    explicit
    BackwardDifferentiationFormula(const unsigned num_steps)
        : _num_steps(num_steps)
    {
        // TODO: assert 0 < num_steps <= 6
        _xs_old.reserve(num_steps);
    }

    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _xs_old.push_back(x);
    }

    void pushState(const double t, Vector const& x, InternalMatrixStorage const&) override
    {
        (void) t;

        // until _xs_old is filled, lower-order BDF formulas are used.
        if (_xs_old.size() < _num_steps) {
            _xs_old.push_back(x);
        } else {
            _xs_old[_offset] = x;
            _offset = (_offset+1) % _num_steps;
        }
    }

    void setCurrentTime(const double t, const double delta_t) override {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override {
        return _t;
    }

    double getCurrentXWeight() const override {
        auto const k = eff_num_steps();
        return detail::BDF_Coeffs[k-1][0] / _delta_t;
    }

    void getWeightedOldX(Vector& y) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const k = eff_num_steps();
        auto const*const BDFk = detail::BDF_Coeffs[k-1];

        // compute linear combination \sum_{i=0}^{k-1} BDFk_{k-i} \cdot x_{n+i}
        BLAS::copy(_xs_old[_offset], y); // _xs_old[offset] = x_n
        BLAS::scale(y, BDFk[k]);

        for (unsigned i=1; i<k; ++i) {
            auto const off = (_offset + i) % k;
            BLAS::axpy(y, BDFk[k-i], _xs_old[off]);
        }

        BLAS::scale(y, 1.0/_delta_t);
    }

private:
    unsigned eff_num_steps() const { return _xs_old.size(); }

    const unsigned _num_steps;
    double _t = 9999.9999;
    double _delta_t = 8888.8888;

    std::vector<Vector> _xs_old;
    unsigned _offset = 0;
};

//! @}

}
