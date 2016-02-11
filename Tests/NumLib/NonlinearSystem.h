#pragma once

#include "ODETypes.h"


template<NonlinearSolverTag NLTag>
class NonlinearSystem;

template<>
class NonlinearSystem<NonlinearSolverTag::Newton>
{
public:
    virtual void assembleResidualNewton(Vector const& x) = 0;
    virtual void assembleJacobian(Vector const& x) = 0;
    virtual Vector getResidual(Vector const& x) = 0;
    virtual Matrix getJacobian() = 0;

    virtual bool isLinear() const = 0;

    virtual ~NonlinearSystem() = default;
};

template<>
class NonlinearSystem<NonlinearSolverTag::Picard>
{
public:
    virtual void assembleMatricesPicard(Vector const& x) = 0;
    virtual void getA(Matrix& A) = 0;
    virtual void getRhs(Vector& rhs) = 0;

    virtual bool isLinear() const = 0;

    virtual ~NonlinearSystem() = default;
};
