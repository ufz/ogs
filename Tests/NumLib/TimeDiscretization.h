#pragma once

#include "ODETypes.h"

class ITimeDiscretization
{
public:
    virtual void setInitialState(const double t0, Vector const& x) = 0;
    virtual void pushState(const double t, Vector const& x) = 0;
    virtual void pushMatrices() = 0;
    virtual void setCurrentTime(const double t, const double delta_t) = 0;
    virtual double getCurrentTime() const = 0; // get time used for assembly

    // \dot x === alpha * x - x_old
    virtual double getCurrentXWeight() = 0; // = alpha
    virtual Vector getWeightedOldX() = 0; // = x_old

    ~ITimeDiscretization() = default;
};


class BackwardEuler : ITimeDiscretization
{
public:
    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _x_old = x;
    }

    void pushState(const double t, Vector const& x) override {
        (void) t;
        _x_old = x;
    }

    void pushMatrices() override {}

    void setCurrentTime(const double t, const double delta_t) override {
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override {
        return _t;
    }

    double getCurrentXWeight() override {
        return 1.0/_delta_t;
    }

    Vector getWeightedOldX() override {
        return _x_old / _delta_t;
    }

private:
    double _t = 9999.9999;
    double _delta_t = 8888.8888;
    Vector _x_old;
};


class ForwardEuler : ITimeDiscretization
{
public:
    void setInitialState(const double t0, Vector const& x) override {
        _t = t0;
        _t_old = t0;
        _x_old = x;
    }

    void pushState(const double t, Vector const& x) override {
        (void) t;
        _x_old = x;
    }

    void pushMatrices() override {}

    void setCurrentTime(const double t, const double delta_t) override {
        _t_old = _t;
        _t = t;
        _delta_t = delta_t;
    }

    double getCurrentTime() const override {
        return _t_old; // forward Euler does assembly at the preceding timestep
    }

    double getCurrentXWeight() override {
        return 1.0/_delta_t;
    }

    Vector getWeightedOldX() override {
        return _x_old / _delta_t;
    }

private:
    double _t = 9999.9999;
    double _t_old = 7777.7777;
    double _delta_t = 8888.8888;
    Vector _x_old;
};
