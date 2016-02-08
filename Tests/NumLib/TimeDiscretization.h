#pragma once

#include "ODETypes.h"

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

    // Forward Euler is linear, other schemes not.
    virtual bool isLinearTimeDisc() const { return false; }

    // for Crank-Nicolson
    virtual bool needsPreload() const { return false; }

    // \dot x === alpha * x - x_old
    virtual double getCurrentXWeight() const = 0; // = alpha
    virtual Vector getWeightedOldX() const = 0; // = x_old

    ~ITimeDiscretization() = default;

protected:

    // for Crank-Nicolson
    virtual void adjustMatrix(Matrix& A) const { (void) A; }

    // for Crank-Nicolson
    virtual void adjustRhs(Vector& rhs) const { (void) rhs; }

    // for Crank-Nicolson
    virtual void adjustResidual(Vector const& x, Vector& res) const {
        (void) res; (void) x;
    }
};


class BackwardEuler : public ITimeDiscretization
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


class ForwardEuler : public ITimeDiscretization
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


class CrankNicolson : public ITimeDiscretization
{
public:
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

protected:
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
    const double _theta = 0.5;
    double _t = 9999.9999;
    double _delta_t = 8888.8888;
    Vector _x_old;

    Matrix _M_bar;
    Vector _b_bar;
};
