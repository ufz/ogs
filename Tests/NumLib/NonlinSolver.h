#pragma once

#include <Eigen/SparseCore>
#include <logog/include/logog.hpp>

#include "ODETypes.h"

enum class NonlinearSolverTag : bool { Picard, Newton };

template<NonlinearSolverTag NLTag>
class INonlinearSystem;

template<>
class INonlinearSystem<NonlinearSolverTag::Newton>
{
public:
    virtual void assembleResidualNewton(Vector const& x) = 0;
    virtual void assembleJacobian(Vector const& x) = 0;
    virtual Vector getResidual(Vector const& x) = 0;
    virtual Matrix getJacobian() = 0;

    virtual bool isLinear() const = 0;

    virtual ~INonlinearSystem() = default;
};

template<>
class INonlinearSystem<NonlinearSolverTag::Picard>
{
public:
    virtual void assembleMatricesPicard(Vector const& x) = 0;
    virtual Matrix getA() = 0;
    virtual Vector getRhs() = 0;

    virtual bool isLinear() const = 0;

    virtual ~INonlinearSystem() = default;
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

    // for Crank-Nicolson
    void assemble(INonlinearSystem<NonlinearSolverTag::Newton>& sys, Vector& x);

    void solve(INonlinearSystem<NonlinearSolverTag::Newton>& sys, Vector& x);

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

    // for Crank-Nicolson
    void assemble(INonlinearSystem<NonlinearSolverTag::Picard>& sys, Vector& x);

    void solve(INonlinearSystem<NonlinearSolverTag::Picard>& sys, Vector& x);

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
                                  const double dx_dx,
                                  Matrix& Jac) = 0;
};
