/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TESTS_NUMLIB_ODES_H
#define TESTS_NUMLIB_ODES_H

#include "MathLib/LinAlg/BLAS.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/ODESolver/ODESystem.h"

// debug
//#include <iostream>


template<typename Matrix, typename Vector,
         template<typename /*Matrix*/, typename /*Vector*/> class Ode>
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
        MathLib::setMatrix(M, { 1.0, 0.0,  0.0, 1.0 });
        MathLib::setMatrix(K, { 0.0, 1.0, -1.0, 0.0 });

        MathLib::setVector(b, { 0.0, 0.0 });
    }

    void assembleJacobian(const double /*t*/, const Vector &/*x*/, Vector const& /*xdot*/,
                          const double dxdot_dx,  const Matrix& M,
                          const double dx_dx, const Matrix& K,
                          Matrix &Jac) override
    {
        namespace BLAS = MathLib::BLAS;

        // compute Jac = M*dxdot_dx + dx_dx*K
        BLAS::copy(M, Jac);
        BLAS::scale(Jac, dxdot_dx);
        if (dx_dx != 0.0)
            BLAS::axpy(Jac, dx_dx, K);
    }

    MathLib::MatrixSpecifications getMatrixSpecifications() const override
    {
        return { N, N, nullptr, nullptr, nullptr };
    }

    bool isLinear() const override
    {
        return true;
    }

    std::size_t const N = 2;
};

template<typename Matrix, typename Vector>
class ODETraits<Matrix, Vector, ODE1>
{
public:
    static void setIC(Vector& x0)
    {
        MathLib::setVector(x0, { 1.0, 0.0 });
        MathLib::BLAS::finalizeAssembly(x0);
    }

    static Vector solution(const double t)
    {
        Vector v(2);
        MathLib::setVector(v, { cos(t), sin(t) });
        MathLib::BLAS::finalizeAssembly(v);
        return v;
    }

    static const double t0;
    static const double t_end;
};

template<typename Matrix, typename Vector>
const double ODETraits<Matrix, Vector, ODE1>::t0 = 0.0;

template<typename Matrix, typename Vector>
const double ODETraits<Matrix, Vector, ODE1>::t_end = 30.0;
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
        MathLib::setMatrix(M, { 1.0 });
        MathLib::setMatrix(K, { x[0] });
        MathLib::setVector(b, { 0.0 });
    }

    void assembleJacobian(const double /*t*/, const Vector &x, Vector const& /*xdot*/,
                          const double dxdot_dx, Matrix const& M,
                          const double dx_dx, Matrix const& K,
                          Matrix &Jac) override
    {
        namespace BLAS = MathLib::BLAS;

        // compute Jac = M*dxdot_dx + dK_dx + dx_dx*K
        BLAS::copy(M, Jac);
        BLAS::scale(Jac, dxdot_dx);

        MathLib::addToMatrix(Jac, { x[0] }); // add dK_dx

        if (dx_dx != 0.0)
        {
            BLAS::finalizeAssembly(Jac);
            BLAS::axpy(Jac, dx_dx, K);
        }
    }

    MathLib::MatrixSpecifications getMatrixSpecifications() const override
    {
        return { N, N, nullptr, nullptr, nullptr };
    }

    bool isLinear() const override
    {
        return false;
    }

    std::size_t const N = 1;
};

template<typename Matrix, typename Vector>
class ODETraits<Matrix, Vector, ODE2>
{
public:
    static void setIC(Vector& x0)
    {
        MathLib::setVector(x0, { 1.0 });
        MathLib::BLAS::finalizeAssembly(x0);
    }

    static Vector solution(const double t)
    {
        Vector v(1);
        MathLib::setVector(v, { 1.0/t });
        MathLib::BLAS::finalizeAssembly(v);
        return v;
    }

    static const double t0;
    static const double t_end;
};

template<typename Matrix, typename Vector>
const double ODETraits<Matrix, Vector, ODE2>::t0 = 1.0;

template<typename Matrix, typename Vector>
const double ODETraits<Matrix, Vector, ODE2>::t_end = 2.0;
// ODE 2 end //////////////////////////////////////////////////////

// ODE 3 //////////////////////////////////////////////////////////
//
// TODO Check that the ODE is really solved by the given solution,
//      i.e., residual is zero.
//      Check (by numerical differentiation that the Jacobian is correct.
//
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

        MathLib::
        setMatrix(M, {       t*y, 1.0,     0.0,
                             0.0,  -t,     t*y,
                       omega*x*t, 0.0, omega*x });
        MathLib::
        setMatrix(K, {             y,   1.0/t,                       -y,
                       omega*omega/y,    -0.5,                      0.0,
                        -0.5*omega*z, y/omega, -(1.0/omega/t+omega)*y*z });
        MathLib::
        setVector(b, { 0.0,
                       0.5/t,
                       0.5*omega*x*z + omega/t });
    }

    void assembleJacobian(const double t, const Vector& x_curr, Vector const& xdot,
                          const double dxdot_dx, Matrix const& M,
                          const double dx_dx, Matrix const& K,
                          Matrix &Jac) override
    {
        auto const x = x_curr[0];
        auto const y = x_curr[1];
        auto const z = x_curr[2];

        auto const dx = xdot[0];
        auto const dz = xdot[2];

        namespace BLAS = MathLib::BLAS;

        // Compute Jac = M dxdot/dx + dM/dx xdot + K dx/dx + dK/dx x - db/dx

        BLAS::copy(M, Jac);
        BLAS::scale(Jac, dxdot_dx); // Jac = M * dxdot_dx

        /* dx_dx == 0 holds if and only if the ForwardEuler scheme is used.
         *
         * In that case M, K, and b are assembled at the preceding timestep.
         * Thus their derivatices w.r.t. x of the current timestep vanishes
         * and does not contribute to the Jacobian.
         *
         * TODO: Maybe if relaxation is applied, dx_dx takes values from the
         *       interval [ 0.0, 1.0 ]
         */
        if (dx_dx != 0.0)
        {
            // in this block it is assumed that dx_dx == 1.0

            // add dM/dx \cdot \dot x
            MathLib::
            addToMatrix(Jac, {                 0.0, t*dx, 0.0,
                                               0.0, t*dz, 0.0,
                               omega*t*dx+omega*dz,  0.0, 0.0 });

            BLAS::axpy(Jac, dx_dx, K); // add K \cdot dx_dx

            // add dK/dx \cdot \dot x
            MathLib::
            addToMatrix(Jac, { 0.0, x-z, 0.0,
                               0.0, -omega*omega/y/y*x, 0.0,
                               0.0, y/omega-(1.0/omega/t+omega)*z*z, // -->
                               /* --> */  -0.5*omega*x - (1.0/omega/t+omega)*y*z });

            // add -db/dx
            MathLib::
            addToMatrix(Jac, {          0.0, 0.0,          0.0,
                                        0.0, 0.0,          0.0,
                               -0.5*omega*z, 0.0, -0.5*omega*z });
        }

        // Eigen::MatrixXd J(Jac.getRawMatrix());

        // INFO("t: %e, x: %e, y: %e, z: %e, dxdot/dx: %e", t, x, y, z, dxdot_dx);
        // std::cout << "Jacobian:\n" << J << "\n";

        // INFO("Det J: %e <<<", J.determinant());
    }

    MathLib::MatrixSpecifications getMatrixSpecifications() const override
    {
        return { N, N, nullptr, nullptr, nullptr };
    }

    bool isLinear() const override
    {
        return false;
    }

    std::size_t const N = 3;
    static const double omega;
};

template<typename Matrix, typename Vector>
const double ODE3<Matrix, Vector>::omega = 15.0;


template<typename Matrix, typename Vector>
class ODETraits<Matrix, Vector, ODE3>
{
public:
    static void setIC(Vector& x0)
    {
        auto const omega = ODE3<Matrix, Vector>::omega;

        MathLib::setVector(x0, { sin(omega*t0)/omega/t0,
                                 1.0/t0,
                                 cos(omega*t0) });
        MathLib::BLAS::finalizeAssembly(x0);

        // std::cout << "IC:\n" << Eigen::VectorXd(x0.getRawVector()) << "\n";
    }

    static Vector solution(const double t)
    {
        auto const omega = ODE3<Matrix, Vector>::omega;

        Vector v(3);
        MathLib::setVector(v, { sin(omega*t)/omega/t,
                                1.0/t,
                                cos(omega*t) });
        MathLib::BLAS::finalizeAssembly(v);
        return v;
    }

    static const double t0;
    static const double t_end;
};

template<typename Matrix, typename Vector>
const double ODETraits<Matrix, Vector, ODE3>::t0 = 1.0;

template<typename Matrix, typename Vector>
const double ODETraits<Matrix, Vector, ODE3>::t_end = 3.0;
// ODE 3 end //////////////////////////////////////////////////////

#endif // TESTS_NUMLIB_ODES_H
