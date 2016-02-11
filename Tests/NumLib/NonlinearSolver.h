#pragma once

#include <logog/include/logog.hpp>

#include "ODETypes.h"
#include "NonlinearSystem.h"


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
    void assemble(NonlinearSystem<NonlinearSolverTag::Newton>& sys, Vector& x);

    void solve(NonlinearSystem<NonlinearSolverTag::Newton>& sys, Vector& x);

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
    void assemble(NonlinearSystem<NonlinearSolverTag::Picard>& sys, Vector& x);

    void solve(NonlinearSystem<NonlinearSolverTag::Picard>& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _x_new;
};
