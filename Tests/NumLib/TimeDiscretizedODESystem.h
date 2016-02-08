#pragma once

#include "NonlinSolver.h"
#include "TimeDiscretization.h"

template<NonlinearSolverTag NLTag, typename TimeDisc>
struct TimeDiscretizedODESystem;

template<typename TimeDisc>
class TimeDiscretizedODESystem<NonlinearSolverTag::Newton, TimeDisc> final
        : public INonlinearSystemNewton
        , public IParabolicEquation
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& ode,
                             TimeDisc& time_discretization)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _Jac(ode.getMatrixSize(), ode.getMatrixSize())
        , _M(_Jac)
        , _K(_Jac)
        , _b(ode.getMatrixSize())
    {}

    /// begin INonlinearSystemNewton

    void assembleResidualNewton(const Vector &x) override
    {
        auto const t = _time_disc.getCurrentTime();
        _ode.assemble(t, x, _M, _K, _b);
    }

    void assembleJacobian(const Vector &x) override
    {
        auto const t = _time_disc.getCurrentTime();
        auto const dxdot_dx = _time_disc.getCurrentXWeight();
        _ode.assembleJacobian(t, x, dxdot_dx, _Jac);
        _time_disc.adjustMatrix(_Jac);
    }

    Vector getResidual(Vector const& x) override
    {
        auto const alpha = _time_disc.getCurrentXWeight();
        auto const x_old = _time_disc.getWeightedOldX();
        Vector res = _M * (alpha*x - x_old) + _K*x - _b;
        _time_disc.adjustResidual(x, res);
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

    TimeDisc& getTimeDiscretization() {
        return _time_disc;
    }

private:
    // from ITimeDiscretization
    virtual void getMatrices(Matrix const*& M, Matrix const*& K,
                             Vector const*& b) const override
    {
        M = &_M;
        K = &_K;
        b = &_b;
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& _ode;
    TimeDisc& _time_disc;

    Matrix _Jac;
    Matrix _M;
    Matrix _K;
    Vector _b;
};

template<typename TimeDisc>
class TimeDiscretizedODESystem<NonlinearSolverTag::Picard, TimeDisc> final
        : public INonlinearSystemPicard
        , public IParabolicEquation
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& ode,
                             TimeDisc& time_discretization)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _M(ode.getMatrixSize(), ode.getMatrixSize())
        , _K(_M)
        , _b(ode.getMatrixSize())
    {}

    /// begin INonlinearSystemNewton

    void assembleMatricesPicard(const Vector &x) override
    {
        auto const t = _time_disc.getCurrentTime();
        _ode.assemble(t, x, _M, _K, _b);
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

    TimeDisc& getTimeDiscretization() {
        return _time_disc;
    }

private:
    // from ITimeDiscretization
    virtual void getMatrices(Matrix const*& M, Matrix const*& K,
                             Vector const*& b) const override
    {
        M = &_M;
        K = &_K;
        b = &_b;
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& _ode;
    TimeDisc& _time_disc;

    Matrix _M;
    Matrix _K;
    Vector _b;
};
