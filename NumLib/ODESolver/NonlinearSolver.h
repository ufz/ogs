#pragma once

#include <logog/include/logog.hpp>

#include "Types.h"
#include "NonlinearSystem.h"


namespace NumLib
{

template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class NonlinearSolver;

template<typename Matrix, typename Vector>
class NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Newton> final
{
public:
    using System = NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton>;

    explicit
    NonlinearSolver(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    // for Crank-Nicolson
    void assemble(System& sys, Vector& x);

    bool solve(System& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _minus_delta_x;
};

template<typename Matrix, typename Vector>
class NonlinearSolver<Matrix, Vector, NonlinearSolverTag::Picard> final
{
public:
    using System = NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard>;

    explicit
    NonlinearSolver(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    // for Crank-Nicolson
    void assemble(System& sys, Vector& x);

    bool solve(System& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _x_new;
};

}

#include "NonlinearSolver-impl.h"
