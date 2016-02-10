#pragma once

#include "ODETypes.h"

#include <vector>

// debugging
#include <iostream>

class IParabolicEquation
{
public:
    // needed for Crank-Nicolson
    virtual void getMatrices(Matrix const*& M, Matrix const*& K,
                             Vector const*& b) const = 0;
};

class ITimeDiscretization
{
public:
    virtual void setInitialState(const double t0, Vector const& x) = 0;
    virtual void pushState(const double t, Vector const& x,
                           IParabolicEquation const& eq) = 0;

    virtual void setCurrentTime(const double t, const double delta_t) = 0;
    virtual double getCurrentTime() const = 0; // get time used for assembly

    // \dot x === alpha * x - x_old
    virtual double getCurrentXWeight() const = 0; // = alpha
    virtual Vector getWeightedOldX() const = 0; // = x_old

    ~ITimeDiscretization() = default;


    // Forward Euler is linear, other schemes not.
    virtual bool isLinearTimeDisc() const { return false; }

    // Forward Euler overrides this.
    // Caution: This is not the x with which you want to compute \dot x
    virtual Vector const& getCurrentX(Vector const& x_at_new_timestep) const
    {
        return x_at_new_timestep;
    }

    // for Crank-Nicolson
    virtual bool needsPreload() const { return false; }

    // for Crank-Nicolson
    virtual void adjustMatrix(Matrix& A) const { (void) A; }

    // for Crank-Nicolson
    virtual void adjustRhs(Vector& rhs) const { (void) rhs; }

    // for Crank-Nicolson
    virtual void adjustResidual(Vector const& x, Vector& res) const {
        (void) res; (void) x;
    }
};


class BackwardEuler final : public ITimeDiscretization
{
public:
    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _x_old = x;
    }

    void pushState(const double t, Vector const& x, IParabolicEquation const&) override
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

    Vector getWeightedOldX() const override {
        return _x_old / _delta_t;
    }

private:
    double _t = 9999.9999;
    double _delta_t = 8888.8888;
    Vector _x_old;
};


class ForwardEuler final : public ITimeDiscretization
{
public:
    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _t_old = t0;
        _x_old = x;
    }

    void pushState(const double t, Vector const& x, IParabolicEquation const&) override
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

    Vector getWeightedOldX() const override {
        return _x_old / _delta_t;
    }

    bool isLinearTimeDisc() const override {
        return true;
    }

private:
    double _t = 9999.9999;
    double _t_old = 7777.7777;
    double _delta_t = 8888.8888;
    Vector _x_old;
};


class CrankNicolson final : public ITimeDiscretization
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

    void pushState(const double t, Vector const& x, IParabolicEquation const& eq) override
    {
        (void) t;
        _x_old = x;

        Matrix const* M;
        Matrix const* K;
        Vector const* b;
        eq.getMatrices(M, K, b);

        _M_bar = (1-_theta) * (*M);
        _b_bar = (1-_theta) * ((*K)*_x_old - (*b));
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

    Vector getWeightedOldX() const override {
        return _x_old / _delta_t;
    }

    bool needsPreload() const override {
        return true;
    }

    void adjustMatrix(Matrix& A) const override
    {
        auto const alpha = getCurrentXWeight();
        A *= _theta;
        A += alpha * _M_bar;
    }

    void adjustRhs(Vector& rhs) const override
    {
        rhs *= _theta;
        rhs += _M_bar * getWeightedOldX() - _b_bar;
    }

    void adjustResidual(Vector const& x, Vector& res) const override
    {
        auto const alpha = getCurrentXWeight();
        Vector xdot = alpha * x - getWeightedOldX();
        res *= _theta;
        res += _M_bar * xdot + _b_bar;
    }

private:
    const double _theta = 555.555;
    double _t = 9999.9999;
    double _delta_t = 8888.8888;
    Vector _x_old;

    Matrix _M_bar;
    Vector _b_bar;
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


class BackwardDifferentiationFormula final : public ITimeDiscretization
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

    void pushState(const double t, Vector const& x, IParabolicEquation const&) override
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

    Vector getWeightedOldX() const override {
        auto const k = eff_num_steps();
        auto const*const BDFk = detail::BDF_Coeffs[k-1];

        // compute linear combination \sum_{i=0}^{k-1} BDFk_{k-i} \cdot x_{n+i}
        Vector y(BDFk[k]*_xs_old[_offset]); // _xs_old[offset] = x_n

        for (unsigned i=1; i<k; ++i) {
            auto const off = (_offset + i) % k;
            y += BDFk[k-i] * _xs_old[off];
        }

        y /= _delta_t;

        return y;
    }

private:
    unsigned eff_num_steps() const { return _xs_old.size(); }

    const unsigned _num_steps;
    double _t = 9999.9999;
    double _delta_t = 8888.8888;

    std::vector<Vector> _xs_old;
    unsigned _offset = 0;
};
