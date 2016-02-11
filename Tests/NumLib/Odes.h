#pragma once

#include "NonlinSolver.h"

template<typename Ode>
class OdeTraits;

// ODE 1 //////////////////////////////////////////////////////////
class Ode1 final
        : public FirstOrderImplicitQuasilinearODESystem<NonlinearSolverTag::Newton>
{
public:
    void assemble(const double t, Vector const& x,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        (void) t; (void) x;

        Eigen::Matrix2d m;
        m << 1.0, 0.0,
             0.0, 1.0;
        M = m.sparseView();

        Eigen::Matrix2d k;
        k << 0.0, 1.0,
            -1.0, 0.0;
        K = k.sparseView();

        Eigen::Vector2d b_;
        b_ << 0.0, 0.0;
        b = b_;
    }

    void assembleJacobian(const double t, const Vector &x,
                          const double dxdot_dx,  const double dx_dx,
                          Matrix &Jac)
    {
        (void) t; (void) x;

        Eigen::Matrix2d m;
        m << 1.0, 0.0,
             0.0, 1.0;

        Eigen::Matrix2d k;
        k << 0.0, 1.0,
            -1.0, 0.0;

        Jac = (m*dxdot_dx + dx_dx*k).sparseView();
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

template<>
class OdeTraits<Ode1>
{
public:
    static void setIC(Vector& x0)
    {
        x0.resize(2);
        x0 << 1.0, 0.0;
    }

    static Vector solution(const double t)
    {
        Eigen::Vector2d v;
        v << cos(t), sin(t);
        return v;
    }

    static constexpr double t0    = 0.0;
    static constexpr double t_end = 30.0;
};
// ODE 1 end //////////////////////////////////////////////////////

// ODE 2 //////////////////////////////////////////////////////////
class Ode2 final
        : public FirstOrderImplicitQuasilinearODESystem<NonlinearSolverTag::Newton>
{
public:
    void assemble(const double t, Vector const& x,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        (void) t;
        M.coeffRef(0, 0) = 1.0;
        K.coeffRef(0, 0) = x[0];
        b.coeffRef(0, 0) = 0.0;
    }

    void assembleJacobian(const double t, const Vector &x,
                          const double dxdot_dx, const double dx_dx,
                          Matrix &Jac)
    {
        (void) t;
        Jac.coeffRef(0, 0) = dxdot_dx + x[0] + x[0]*dx_dx;
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

template<>
class OdeTraits<Ode2>
{
public:
    static void setIC(Vector& x0)
    {
        x0.resize(1);
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
