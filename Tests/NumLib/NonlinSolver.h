#pragma once

#include <Eigen/SparseCore>
#include <logog/include/logog.hpp>

#include "TimeDiscretization.h"


enum class NonlinearSolverTag : bool { Picard, Newton };


class INonlinearSystem
{
public:
    virtual bool isLinear() const = 0;

    virtual ~INonlinearSystem() = default;
};

class INonlinearSystemNewton : public INonlinearSystem
{
public:
    virtual void assembleResidualNewton(Vector const& x) = 0;
    virtual void assembleJacobian(Vector const& x) = 0;
    virtual Vector getResidual(Vector const& x) = 0;
    virtual Matrix getJacobian() = 0;

    virtual ~INonlinearSystemNewton() = default;
};

class INonlinearSystemPicard : public INonlinearSystem
{
public:
    virtual void assembleMatricesPicard(Vector const& x) = 0;
    virtual Matrix getA() = 0;
    virtual Vector getRhs() = 0;

    virtual ~INonlinearSystemPicard() = default;
};


template<NonlinearSolverTag NLTag>
class NonlinearSolver;

template<>
class NonlinearSolver<NonlinearSolverTag::Newton> final
{
public:
    explicit
    NonlinearSolver(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    void solve(INonlinearSystemNewton& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _minus_delta_x;
};

template<>
class NonlinearSolver<NonlinearSolverTag::Picard> final
{
public:
    explicit
    NonlinearSolver(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    void solve(INonlinearSystemPicard& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _x_new;
};


template<NonlinearSolverTag NLTag>
class IFirstOrderImplicitOde;

template<>
class IFirstOrderImplicitOde<NonlinearSolverTag::Picard>
{
public:
    virtual bool isLinear() const = 0;
    virtual IndexType getMatrixSize() const = 0;

    virtual void assemble(const double t, Vector const& x,
                          Matrix& M, Matrix& K, Vector& b) = 0;


    virtual ~IFirstOrderImplicitOde() = default;
};

template<>
class IFirstOrderImplicitOde<NonlinearSolverTag::Newton>
        : public IFirstOrderImplicitOde<NonlinearSolverTag::Picard>
{
public:
    virtual void assembleJacobian(const double t, Vector const& x,
                                  const double dxdot_dx,
                                  Matrix& Jac) = 0;
};




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
    }

    Vector getResidual(Vector const& x) override
    {
        auto const alpha = TimeDisc::getCurrentXWeight();
        auto const x_old = TimeDisc::getWeightedOldX();
        return _M * (alpha*x - x_old) + _K*x - _b;
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
        return _M * TimeDisc::getCurrentXWeight() + _K;
    }

    Vector getRhs() override
    {
        return _b + _M * TimeDisc::getWeightedOldX();
    }

    bool isLinear() const override
    {
        return TimeDisc::isLinearTimeDisc() || _ode.isLinear();
    }

private:
    IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& _ode;

    Matrix _M;
    Matrix _K;
    Vector _b;
};
