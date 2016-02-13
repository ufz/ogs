#pragma once

#include "NumLib/ODESolver/ODESystem.h"

// debug
#include <iostream>


template<typename Matrix, typename Vector,
         template<typename /*Matrix*/, typename /*Vector*/> typename Ode>
class ODETraits;

// ODE 1 //////////////////////////////////////////////////////////
template<typename Matrix, typename Vector>
class ODE1 final
        : public NumLib::ODESystem<Matrix, Vector,
                           NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
                           NumLib::NonlinearSolverTag::Newton>
{
public:
    void assemble(const double /*t*/, Vector const& /*x*/,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        setMatrix(M, N, N, { 1.0, 0.0,  0.0, 1.0 });
        setMatrix(K, N, N, { 0.0, 1.0, -1.0, 0.0 });

        setVector(b, { 0.0, 0.0 });
    }

    void assembleJacobian(const double /*t*/, const Vector &/*x*/, Vector const& /*xdot*/,
                          const double dxdot_dx,  const double dx_dx,
                          Matrix &Jac)
    {
        Eigen::MatrixXd m(N, N);
        m << 1.0, 0.0,
             0.0, 1.0;

        Eigen::MatrixXd k(N, N);
        k << 0.0, 1.0,
            -1.0, 0.0;

        setMatrix(Jac, m*dxdot_dx + dx_dx*k);
    }

    NumLib::IndexType getMatrixSize() const override
    {
        return N;
    }

    bool isLinear() const override
    {
        return true;
    }

    const unsigned N = 2;
};

template<typename Matrix, typename Vector>
class ODETraits<Matrix, Vector, ODE1>
{
public:
    static void setIC(Vector& x0)
    {
        setVector(x0, { 1.0, 0.0 });
    }

    static Vector solution(const double t)
    {
        Vector v(2);
        setVector(v, { cos(t), sin(t) });
        return v;
    }

    static constexpr double t0    = 0.0;
    static constexpr double t_end = 30.0;
};
// ODE 1 end //////////////////////////////////////////////////////

// ODE 2 //////////////////////////////////////////////////////////
template<typename Matrix, typename Vector>
class ODE2 final
        : public NumLib::ODESystem<Matrix, Vector,
                           NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
                           NumLib::NonlinearSolverTag::Newton>
{
public:
    void assemble(const double /*t*/, Vector const& x,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        setMatrix(M, N, N, { 1.0 });
        setMatrix(K, N, N, { x[0] });
        setVector(b, { 0.0 });
    }

    void assembleJacobian(const double /*t*/, const Vector &x, Vector const& /*xdot*/,
                          const double dxdot_dx, const double dx_dx,
                          Matrix &Jac)
    {
        setMatrix(Jac, N, N, { dxdot_dx + x[0] + x[0]*dx_dx });
    }

    NumLib::IndexType getMatrixSize() const override
    {
        return N;
    }

    bool isLinear() const override
    {
        return false;
    }

    const unsigned N = 1;
};

template<typename Matrix, typename Vector>
class ODETraits<Matrix, Vector, ODE2>
{
public:
    static void setIC(Vector& x0)
    {
        setVector(x0, { 1.0 });
    }

    static Vector solution(const double t)
    {
        Vector v(1);
        setVector(v, { 1.0/t });
        return v;
    }

    static constexpr double t0    = 1.0;
    static constexpr double t_end = 2.0;
};
// ODE 2 end //////////////////////////////////////////////////////

// ODE 3 //////////////////////////////////////////////////////////
template<typename Matrix, typename Vector>
class ODE3 final
        : public NumLib::ODESystem<Matrix, Vector,
                           NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
                           NumLib::NonlinearSolverTag::Newton>
{
public:
    void assemble(const double t, Vector const& x_curr,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        auto const x = x_curr[0];
        auto const y = x_curr[1];
        auto const z = x_curr[2];

        setMatrix(M, N, N, {       t*y, 1.0,     0.0,
                                   0.0,  -t,     t*y,
                             omega*x*t, 0.0, omega*x });

        setMatrix(K, N, N, {             y,   1.0/t,                       -y,
                             omega*omega/y,    -0.5,                      0.0,
                              -0.5*omega*z, y/omega, -(1.0/omega/t+omega)*y*z });

        b[0] = 0.0;
        b[1] = 0.5/t;
        b[2] = 0.5*omega*x*z + omega/t;
    }

    void assembleJacobian(const double t, const Vector& x_curr, Vector const& xdot,
                          const double dxdot_dx,  const double dx_dx,
                          Matrix &Jac)
    {
        auto const x = x_curr[0];
        auto const y = x_curr[1];
        auto const z = x_curr[2];

        auto const dx = xdot[0];
        auto const dz = xdot[1];

        auto const a = dxdot_dx;

        // set J to M \cdot d\dot x/dx
        setMatrix(Jac, N, N, { a*      t*y, a*1.0, a*    0.0,
                               a*      0.0, a* -t, a*    t*y,
                               a*omega*x*t, a*0.0, a*omega*x });

        if (dx_dx != 0.0)
        {
            // add dM/dx \cdot \dot x
            addToMatrix(Jac, N, N, {                 0.0, t*dx, 0.0,
                                                     0.0, t*dz, 0.0,
                                     omega*t*dx+omega*dz,  0.0, 0.0 });

            // add K \cdot dx/dx
            addToMatrix(Jac, N, N, {             y,   1.0/t,                       -y,
                                     omega*omega/y,    -0.5,                      0.0,
                                      -0.5*omega*z, y/omega, -(1.0/omega/t+omega)*y*z });

            // add dK/dx \cdot \dot x
            addToMatrix(Jac, N, N, { 0.0, x-z, 0.0,
                                     0.0, -omega*omega/y/y*x, 0.0,
                                     0.0, y/omega-(1.0/omega/t+omega)*z*z, // -->
                                     /* --> */  -0.5*omega*x - (1.0/omega/t+omega)*y*z });

            // add db/dx
            addToMatrix(Jac, N, N, {          0.0, 0.0,          0.0,
                                              0.0, 0.0,          0.0,
                                     -0.5*omega*z, 0.0, -0.5*omega*z });
        }


        // Eigen::MatrixXd J(Jac.getRawMatrix());

        // INFO("t: %e, x: %e, y: %e, z: %e, dxdot/dx: %e", t, x, y, z, dxdot_dx);
        // std::cout << "Jacobian:\n" << J << "\n";

        // INFO("Det J: %e <<<", J.determinant());
    }

    NumLib::IndexType getMatrixSize() const override
    {
        return N;
    }

    bool isLinear() const override
    {
        return false;
    }

    const unsigned N = 3;
    static constexpr double omega = 15.0;
};

template<typename Matrix, typename Vector>
class ODETraits<Matrix, Vector, ODE3>
{
public:
    static void setIC(Vector& x0)
    {
        auto const omega = ODE3<Matrix, Vector>::omega;

        x0[0] = sin(omega*t0)/omega/t0;
        x0[1] = 1.0/t0;
        x0[2] = cos(omega*t0);

        // std::cout << "IC:\n" << Eigen::VectorXd(x0.getRawVector()) << "\n";
    }

    static Vector solution(const double t)
    {
        auto const omega = ODE3<Matrix, Vector>::omega;

        Vector v(3);
        v[0] = sin(omega*t)/omega/t;
        v[1] = 1.0/t;
        v[2] = cos(omega*t);
        return v;
    }

    static constexpr double t0    = 1.0;
    static constexpr double t_end = 3.0;
};
// ODE 3 end //////////////////////////////////////////////////////
