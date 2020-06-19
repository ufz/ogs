/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/math/constants/constants.hpp>
#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/ODESolver/ODESystem.h"

// debug
//#include <iostream>

template <class Ode>
class ODETraits;

// ODE 1 //////////////////////////////////////////////////////////
class ODE1 final : public NumLib::ODESystem<
                       NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
                       NumLib::NonlinearSolverTag::Newton>
{
public:
    void preAssemble(const double /*t*/, double const /*dt*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void setMKbValues(GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
    {
        MathLib::setMatrix(M, {1.0, 0.0, 0.0, 1.0});
        MathLib::setMatrix(K, {0.0, 1.0, -1.0, 0.0});

        MathLib::setVector(b, {0.0, 0.0});
    }

    void assemble(const double /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& /*x*/,
                  std::vector<GlobalVector*> const& /*xdot*/,
                  int const /*process_id*/, GlobalMatrix& M, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        setMKbValues(M, K, b);
    }

    void assembleWithJacobian(const double /*t*/, double const /*dt*/,
                              std::vector<GlobalVector*> const& /*x_curr*/,
                              std::vector<GlobalVector*> const& /*xdot*/,
                              const double dxdot_dx, const double dx_dx,
                              int const /*process_id*/, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        namespace LinAlg = MathLib::LinAlg;

        setMKbValues(M, K, b);

        // compute Jac = M*dxdot_dx + dx_dx*K
        LinAlg::finalizeAssembly(M);
        LinAlg::copy(M, Jac);
        LinAlg::scale(Jac, dxdot_dx);
        if (dx_dx != 0.0)
        {
            LinAlg::finalizeAssembly(K);
            LinAlg::axpy(Jac, dx_dx, K);
        }
    }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int /*process_id*/) const override
    {
        return { N, N, nullptr, nullptr };
    }

    bool isLinear() const override
    {
        return true;
    }

    std::size_t const N = 2;
};

template <>
class ODETraits<ODE1>
{
public:
    static void setIC(GlobalVector& x0)
    {
        MathLib::setVector(x0, { 1.0, 0.0 });
        MathLib::LinAlg::finalizeAssembly(x0);
    }

    static GlobalVector solution(const double t)
    {
        GlobalVector v(2);
        MathLib::setVector(v, {cos(t), sin(t)});
        MathLib::LinAlg::finalizeAssembly(v);
        return v;
    }

    static const double t0;
    static const double t_end;
};

const double ODETraits<ODE1>::t0 = 0.0;

const double ODETraits<ODE1>::t_end =
    2.0 * 3.141592653589793238462643383279502884197169399375105820974944;
// ODE 1 end //////////////////////////////////////////////////////

// ODE 2 //////////////////////////////////////////////////////////
class ODE2 final : public NumLib::ODESystem<
                       NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
                       NumLib::NonlinearSolverTag::Newton>
{
public:
    void preAssemble(const double /*t*/, double const /*dt*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void setMKbValues(GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
                      GlobalVector& b)
    {
        MathLib::setMatrix(M, {1.0});
        MathLib::LinAlg::setLocalAccessibleVector(x);
        MathLib::setMatrix(K, {x[0]});
        MathLib::setVector(b, {0.0});
    }
    void assemble(const double /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& x,
                  std::vector<GlobalVector*> const& /*xdot*/,
                  int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x[process_id]);
        setMKbValues(*x[process_id], M, K, b);
    }

    void assembleWithJacobian(const double /*t*/, double const /*dt*/,
                              std::vector<GlobalVector*> const& x,
                              std::vector<GlobalVector*> const& /*xdot*/,
                              const double dxdot_dx, const double dx_dx,
                              int const process_id, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x[process_id]);
        setMKbValues(*x[process_id], M, K, b);

        namespace LinAlg = MathLib::LinAlg;

        LinAlg::finalizeAssembly(M);
        // compute Jac = M*dxdot_dx + dK_dx + dx_dx*K
        LinAlg::copy(M, Jac);
        LinAlg::scale(Jac, dxdot_dx);

        MathLib::addToMatrix(Jac, {(*x[process_id])[0]});  // add dK_dx

        if (dx_dx != 0.0)
        {
            LinAlg::finalizeAssembly(K);
            LinAlg::finalizeAssembly(Jac);
            LinAlg::axpy(Jac, dx_dx, K);
        }
    }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int /*process_id*/) const override
    {
        return { N, N, nullptr, nullptr };
    }

    bool isLinear() const override
    {
        return false;
    }

    std::size_t const N = 1;
};

template <>
class ODETraits<ODE2>
{
public:
    static void setIC(GlobalVector& x0)
    {
        MathLib::setVector(x0, { 1.0 });
        MathLib::LinAlg::finalizeAssembly(x0);
    }

    static GlobalVector solution(const double t)
    {
        GlobalVector v(1);
        MathLib::setVector(v, { 1.0/t });
        MathLib::LinAlg::finalizeAssembly(v);
        return v;
    }

    static const double t0;
    static const double t_end;
};

const double ODETraits<ODE2>::t0 = 1.0;

const double ODETraits<ODE2>::t_end = 2.0;
// ODE 2 end //////////////////////////////////////////////////////

// ODE 3 //////////////////////////////////////////////////////////
//
// TODO Check that the ODE is really solved by the given solution,
//      i.e., residual is zero.
//      Check (by numerical differentiation that the Jacobian is correct.
//
class ODE3 final : public NumLib::ODESystem<
                       NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
                       NumLib::NonlinearSolverTag::Newton>
{
public:
    void preAssemble(const double /*t*/, double const /*dt*/,
                     GlobalVector const& /*x*/) override
    {
    }

    void setMKbValues(GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
                      GlobalVector& b)
    {
        auto const u = x[0];
        auto const v = x[1];

        MathLib::setMatrix(M, {u, 2.0 - u, 2.0 - v, v});
        MathLib::setMatrix(K, {2.0 - u - v, -u, v, 2.0 * v - 5.0});
        MathLib::setVector(
            b, {-2 * u + 2.0 * u * v, -2.0 * v * v + 5.0 * v - 4.0});
    }

    void assemble(const double /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& x_curr,
                  std::vector<GlobalVector*> const& /*xdot*/,
                  int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x_curr[process_id]);
        setMKbValues(*x_curr[process_id], M, K, b);
    }

    void assembleWithJacobian(const double /*t*/, double const /*dt*/,
                              std::vector<GlobalVector*> const& x_curr,
                              std::vector<GlobalVector*> const& xdot,
                              const double dxdot_dx, const double dx_dx,
                              int const process_id, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x_curr[process_id]);
        MathLib::LinAlg::setLocalAccessibleVector(*xdot[process_id]);

        setMKbValues(*x_curr[process_id], M, K, b);

        auto const u = (*x_curr[process_id])[0];
        auto const v = (*x_curr[process_id])[1];

        auto const du = (*xdot[process_id])[0];
        auto const dv = (*xdot[process_id])[1];

        namespace LinAlg = MathLib::LinAlg;

        // Compute Jac = M dxdot/dx + dM/dx xdot + K dx/dx + dK/dx x - db/dx

        LinAlg::finalizeAssembly(M);
        LinAlg::copy(M, Jac);
        LinAlg::scale(Jac, dxdot_dx); // Jac = M * dxdot_dx

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
            MathLib::addToMatrix(Jac, {du - dv, 0.0, 0.0, dv - du});

            LinAlg::finalizeAssembly(K);
            LinAlg::finalizeAssembly(Jac);
            LinAlg::axpy(Jac, dx_dx, K); // add K \cdot dx_dx

            // add dK/dx \cdot \dot x
            MathLib::addToMatrix(Jac, {-du - dv, -du, 0.0, du + 2.0 * dv});

            // add -db/dx
            MathLib::addToMatrix(Jac,
                                 {2.0 - 2.0 * v, -2.0 * u, 0.0, 2.0 * v - 5.0});
        }

        // Eigen::MatrixXd J(Jac.getRawMatrix());

        // INFO("t: {:e}, x: {:e}, y: {:e}, z: {:e}, dxdot/dx: {:e}", t, x, y,
        // z, dxdot_dx); std::cout << "Jacobian:\n" << J << "\n";

        // INFO("Det J: {:e} <<<", J.determinant());
    }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int /*process_id*/) const override
    {
        return { N, N, nullptr, nullptr };
    }

    bool isLinear() const override
    {
        return false;
    }

    std::size_t const N = 2;
};

template <>
class ODETraits<ODE3>
{
public:
    static void setIC(GlobalVector& x0)
    {
        MathLib::setVector(x0,
                           {std::sin(2.0 * t0), std::cos(t0) * std::cos(t0)});
        MathLib::LinAlg::finalizeAssembly(x0);

        // std::cout << "IC:\n" << Eigen::VectorXd(x0.getRawVector()) << "\n";
    }

    static GlobalVector solution(const double t)
    {
        GlobalVector v(2);
        MathLib::setVector(v, {std::sin(2.0 * t), std::cos(t) * std::cos(t)});
        MathLib::LinAlg::finalizeAssembly(v);
        return v;
    }

    static const double t0;
    static const double t_end;
};

const double ODETraits<ODE3>::t0 = 0.0;

const double ODETraits<ODE3>::t_end =
    0.5 * boost::math::constants::pi<double>();
// ODE 3 end //////////////////////////////////////////////////////
