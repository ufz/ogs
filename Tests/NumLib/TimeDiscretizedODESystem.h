#pragma once

#include "NonlinSolver.h"

template<NonlinearSolverTag NLTag, typename TimeDisc>
struct TimeDiscretizedODESystem;

template<typename TimeDisc>
class TimeDiscretizedODESystem<NonlinearSolverTag::Newton, TimeDisc> final
        : public TimeDisc
        , public INonlinearSystemNewton
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& ode)
        : _ode(ode)
        , _Jac(ode.getMatrixSize(), ode.getMatrixSize())
        , _M(_Jac)
        , _K(_Jac)
        , _b(ode.getMatrixSize())
    {}

    void assembleResidualNewton(const Vector &x) override
    {
        auto const t = TimeDisc::getCurrentTime();
        _ode.assemble(t, x, _M, _K, _b);
    }

    void assembleJacobian(const Vector &x) override
    {
        auto const t = TimeDisc::getCurrentTime();
        auto const dxdot_dx = TimeDisc::getCurrentXWeight();
        _ode.assembleJacobian(t, x, dxdot_dx, _Jac);
        TimeDisc::adjustMatrix(_Jac);
    }

    Vector getResidual(Vector const& x) override
    {
        auto const alpha = TimeDisc::getCurrentXWeight();
        auto const x_old = TimeDisc::getWeightedOldX();
        Vector res = _M * (alpha*x - x_old) + _K*x - _b;
        TimeDisc::adjustResidual(x, res);
        return res;
    }

    Matrix getJacobian() override
    {
        return _Jac;
    }

    bool isLinear() const override
    {
        return TimeDisc::isLinearTimeDisc() || _ode.isLinear();
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

    Matrix _Jac;
    Matrix _M;
    Matrix _K;
    Vector _b;
};

template<typename TimeDisc>
class TimeDiscretizedODESystem<NonlinearSolverTag::Picard, TimeDisc> final
        : public TimeDisc
        , public INonlinearSystemPicard
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& ode)
        : _ode(ode)
        , _M(ode.getMatrixSize(), ode.getMatrixSize())
        , _K(_M)
        , _b(ode.getMatrixSize())
    {}

    void assembleMatricesPicard(const Vector &x) override
    {
        auto const t = TimeDisc::getCurrentTime();
        _ode.assemble(t, x, _M, _K, _b);
    }

    Matrix getA() override
    {
        Matrix A = _M * TimeDisc::getCurrentXWeight() + _K;
        TimeDisc::adjustMatrix(A);
        return A;
    }

    Vector getRhs() override
    {
        Vector rhs = _b + _M * TimeDisc::getWeightedOldX();
        TimeDisc::adjustRhs(rhs);
        return rhs;
    }

    bool isLinear() const override
    {
        return TimeDisc::isLinearTimeDisc() || _ode.isLinear();
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

    Matrix _M;
    Matrix _K;
    Vector _b;
};
