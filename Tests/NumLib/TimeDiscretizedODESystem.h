#pragma once

#include "NonlinSolver.h"
#include "TimeDiscretization.h"

template<NonlinearSolverTag NLTag>
struct TimeDiscretizedODESystem;

template<>
class TimeDiscretizedODESystem<NonlinearSolverTag::Newton> final
        : public INonlinearSystemNewton
        , public IParabolicEquation
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& ode,
                             ITimeDiscretization& time_discretization)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _Jac(ode.getMatrixSize(), ode.getMatrixSize())
        , _M(_Jac)
        , _K(_Jac)
        , _b(ode.getMatrixSize())
    {}

    /// begin INonlinearSystemNewton

    void assembleResidualNewton(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    void assembleJacobian(const Vector &x_new_timestep) override
    {
        auto const  t        = _time_disc.getCurrentTime();
        auto const& x_curr   = _time_disc.getCurrentX(x_new_timestep);
        auto const  dxdot_dx = _time_disc.getCurrentXWeight();

        _ode.assembleJacobian(t, x_curr, dxdot_dx, _Jac);
        _time_disc.adjustMatrix(_Jac);
    }

    Vector getResidual(Vector const& x_new_timestep) override
    {
        auto const  alpha  = _time_disc.getCurrentXWeight();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);
        auto const  x_old  = _time_disc.getWeightedOldX();
        auto const  x_dot  = alpha*x_new_timestep - x_old;

        Vector res = _M * x_dot + _K*x_curr - _b;
        _time_disc.adjustResidual(x_new_timestep, res);

        return res;
    }

    Matrix getJacobian() override
    {
        return _Jac;
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    /// end INonlinearSystemNewton

    ITimeDiscretization& getTimeDiscretization() {
        return _time_disc;
    }

private:
    // from IParabolicEquation
    virtual void getMatrices(Matrix const*& M, Matrix const*& K,
                             Vector const*& b) const override
    {
        M = &_M;
        K = &_K;
        b = &_b;
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& _ode;
    ITimeDiscretization& _time_disc;

    Matrix _Jac;
    Matrix _M;
    Matrix _K;
    Vector _b;
};

template<>
class TimeDiscretizedODESystem<NonlinearSolverTag::Picard> final
        : public INonlinearSystemPicard
        , public IParabolicEquation
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& ode,
                             ITimeDiscretization& time_discretization)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _M(ode.getMatrixSize(), ode.getMatrixSize())
        , _K(_M)
        , _b(ode.getMatrixSize())
    {}

    /// begin INonlinearSystemNewton

    void assembleMatricesPicard(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    Matrix getA() override
    {
        Matrix A = _M * _time_disc.getCurrentXWeight() + _K;
        _time_disc.adjustMatrix(A);
        return A;
    }

    Vector getRhs() override
    {
        Vector rhs = _b + _M * _time_disc.getWeightedOldX();
        _time_disc.adjustRhs(rhs);
        return rhs;
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    /// end INonlinearSystemPicard

    ITimeDiscretization& getTimeDiscretization() {
        return _time_disc;
    }

private:
    // from IParabolicEquation
    virtual void getMatrices(Matrix const*& M, Matrix const*& K,
                             Vector const*& b) const override
    {
        M = &_M;
        K = &_K;
        b = &_b;
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& _ode;
    ITimeDiscretization& _time_disc;

    Matrix _M;
    Matrix _K;
    Vector _b;
};
