#pragma once

#include "Types.h"


namespace NumLib
{

template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class NonlinearSystem;

template<typename Matrix, typename Vector>
class NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Newton>
{
public:
    virtual void assembleResidualNewton(Vector const& x) = 0;
    virtual void assembleJacobian(Vector const& x) = 0;
    virtual void getResidual(Vector const& x, Vector& res) = 0;
    virtual void getJacobian(Matrix& Jac) = 0;

    virtual bool isLinear() const = 0;

    virtual ~NonlinearSystem() = default;
};

template<typename Matrix, typename Vector>
class NonlinearSystem<Matrix, Vector, NonlinearSolverTag::Picard>
{
public:
    virtual void assembleMatricesPicard(Vector const& x) = 0;
    virtual void getA(Matrix& A) = 0;
    virtual void getRhs(Vector& rhs) = 0;

    virtual bool isLinear() const = 0;

    virtual ~NonlinearSystem() = default;
};

}
