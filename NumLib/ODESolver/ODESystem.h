#pragma once

#include "Types.h"


namespace NumLib
{

//! \addtogroup ODESolver
//! @{

template<typename Matrix, typename Vector, ODESystemTag ODETag, NonlinearSolverTag NLTag>
class ODESystem;

template<typename Matrix, typename Vector>
class ODESystem<Matrix, Vector,
                ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Picard>
{
public:
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;


    virtual bool isLinear() const = 0;
    virtual IndexType getMatrixSize() const = 0;

    virtual void assemble(const double t, Vector const& x,
                          Matrix& M, Matrix& K, Vector& b) = 0;

    virtual ~ODESystem() = default;
};

template<typename Matrix, typename Vector>
class ODESystem<Matrix, Vector,
                ODESystemTag::FirstOrderImplicitQuasilinear,
                NonlinearSolverTag::Newton>
        : public ODESystem<Matrix, Vector,
                           ODESystemTag::FirstOrderImplicitQuasilinear,
                           NonlinearSolverTag::Picard>
{
public:
    virtual void assembleJacobian(const double t, Vector const& x, Vector const& xdot,
                                  const double dxdot_dx, const double dx_dx,
                                  Matrix& Jac) = 0;
};

//! @}

}
