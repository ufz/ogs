#pragma once

#include "ODESystem.h"


template<typename Vector,
         template<typename /*Matrix*/, typename /*Vector*/> typename Ode>
class OdeTraits;

// ODE 1 //////////////////////////////////////////////////////////
template<typename Matrix, typename Vector>
class Ode1 final
        : public ODESystem<Matrix, Vector,
                           ODESystemTag::FirstOrderImplicitQuasilinear,
                           NonlinearSolverTag::Newton>
{
public:
    void assemble(const double /*t*/, Vector const& /*x*/,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        setMatrix(M, 2, 2, { 1.0, 0.0,  0.0, 1.0 });
        setMatrix(K, 2, 2, { 0.0, 1.0, -1.0, 0.0 });

        b[0] = 0.0; b[1] = 0.0;
    }

    void assembleJacobian(const double /*t*/, const Vector &/*x*/,
                          const double dxdot_dx,  const double dx_dx,
                          Matrix &Jac)
    {
        Eigen::Matrix2d m;
        m << 1.0, 0.0,
             0.0, 1.0;

        Eigen::Matrix2d k;
        k << 0.0, 1.0,
            -1.0, 0.0;

        setMatrix(Jac, m*dxdot_dx + dx_dx*k);
    }

    IndexType getMatrixSize() const override
    {
        return 2;
    }

    bool isLinear() const override
    {
        return true;
    }
};

template<typename Vector>
class OdeTraits<Vector, Ode1>
{
public:
    static void setIC(Vector& x0)
    {
        x0[0] = 1.0; x0[1] = 0.0;
    }

    static Vector solution(const double t)
    {
        Vector v(2);
        v[0] = cos(t);
        v[1] = sin(t);
        return v;
    }

    static constexpr double t0    = 0.0;
    static constexpr double t_end = 30.0;
};
// ODE 1 end //////////////////////////////////////////////////////

// ODE 2 //////////////////////////////////////////////////////////
template<typename Matrix, typename Vector>
class Ode2 final
        : public ODESystem<Matrix, Vector,
                           ODESystemTag::FirstOrderImplicitQuasilinear,
                           NonlinearSolverTag::Newton>
{
public:
    void assemble(const double /*t*/, Vector const& x,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        setMatrix(M, 1, 1, { 1.0 });
        setMatrix(K, 1, 1, { x[0] });
        b[0] = 0.0;
    }

    void assembleJacobian(const double /*t*/, const Vector &x,
                          const double dxdot_dx, const double dx_dx,
                          Matrix &Jac)
    {
        setMatrix(Jac, 1, 1, { dxdot_dx + x[0] + x[0]*dx_dx });
    }

    IndexType getMatrixSize() const override
    {
        return 1;
    }

    bool isLinear() const override
    {
        return false;
    }
};

template<typename Vector>
class OdeTraits<Vector, Ode2>
{
public:
    static void setIC(Vector& x0)
    {
        x0[0] = 1.0;
    }

    static Vector solution(const double t)
    {
        Vector v(1);
        v[0] = 1.0 / t;
        return v;
    }

    static constexpr double t0    = 1.0;
    static constexpr double t_end = 2.0;
};
// ODE 2 end //////////////////////////////////////////////////////
