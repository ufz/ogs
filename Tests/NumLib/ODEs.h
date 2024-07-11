/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <numbers>

#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/ODESolver/ODESystem.h"

template <class Ode>
class ODETraits;

inline GlobalVector computeResidumFromMKb(double const dt,
                                          GlobalVector const& x_curr,
                                          GlobalVector const& x_prev,
                                          GlobalMatrix const& M,
                                          GlobalMatrix const& K,
                                          GlobalVector const& b)
{
    using namespace MathLib::LinAlg;
    GlobalVector res(b);

    // res = -M * x_dot - K * x_curr + b
    GlobalVector x_dot;
    copy(x_curr, x_dot);              // x_dot = x
    axpy(x_dot, -1., x_prev);         // x_dot = x - x_prev
    scale(x_dot, 1. / dt);            // x_dot = (x - x_prev)/dt
    matMult(M, x_dot, res);           // res = M*x_dot
    matMultAdd(K, x_curr, res, res);  // res = M*x_dot + K*x
    axpy(res, -1., b);                // res = M*x_dot + K*x - b
    scale(res, -1.);                  // res = -M*x_dot - K*x + b
    return res;
}

inline GlobalMatrix computeJacobianFromMK(double const dt,
                                          GlobalMatrix const& M,
                                          GlobalMatrix const& K)
{
    using namespace MathLib::LinAlg;
    GlobalMatrix J(M);
    // compute J = M*1/dt + K
    copy(M, J);
    scale(J, 1. / dt);
    axpy(J, 1., K);
    finalizeAssembly(J);

    return J;
}

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

    void assembleWithJacobian(const double /*t*/, double const dt,
                              std::vector<GlobalVector*> const& x_curr,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x_curr[process_id]);
        GlobalMatrix M(Jac);
        GlobalMatrix K(Jac);
        setMKbValues(M, K, b);

        MathLib::LinAlg::finalizeAssembly(M);
        MathLib::LinAlg::finalizeAssembly(K);
        MathLib::LinAlg::copy(computeJacobianFromMK(dt, M, K), Jac);
        MathLib::LinAlg::copy(
            computeResidumFromMKb(dt, *x_curr[process_id], *x_prev[process_id],
                                  M, K, b),
            b);
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

    bool requiresNormalization() const override { return false; }

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

const double ODETraits<ODE1>::t_end = 2. * std::numbers::pi;
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

    void assembleWithJacobian(const double /*t*/, double const dt,
                              std::vector<GlobalVector*> const& x_curr,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x_curr[process_id]);
        GlobalMatrix M(Jac);
        GlobalMatrix K(Jac);
        setMKbValues(*x_curr[process_id], M, K, b);

        MathLib::LinAlg::finalizeAssembly(M);
        MathLib::LinAlg::finalizeAssembly(K);
        MathLib::LinAlg::copy(computeJacobianFromMK(dt, M, K), Jac);
        MathLib::LinAlg::copy(
            computeResidumFromMKb(dt, *x_curr[process_id], *x_prev[process_id],
                                  M, K, b),
            b);
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

    bool requiresNormalization() const override { return false; }

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

    void assembleWithJacobian(const double /*t*/, double const dt,
                              std::vector<GlobalVector*> const& x_curr,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalVector& b,
                              GlobalMatrix& Jac) override
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x_curr[process_id]);
        MathLib::LinAlg::setLocalAccessibleVector(*x_prev[process_id]);

        GlobalMatrix M(Jac);
        GlobalMatrix K(Jac);
        setMKbValues(*x_curr[process_id], M, K, b);

        auto const u = (*x_curr[process_id])[0];
        auto const v = (*x_curr[process_id])[1];

        auto const du = (*x_prev[process_id])[0];
        auto const dv = (*x_prev[process_id])[1];

        namespace LinAlg = MathLib::LinAlg;

        // Compute Jac = M dxdot/dx + dM/dx xdot + K dx/dx + dK/dx x - db/dx
        // with dxdot/dx = 1/dt and dx/dx = 1

        LinAlg::finalizeAssembly(M);
        LinAlg::finalizeAssembly(K);
        LinAlg::copy(computeJacobianFromMK(dt, M, K), Jac);
        LinAlg::copy(computeResidumFromMKb(dt, *x_curr[process_id],
                                           *x_prev[process_id], M, K, b),
                     b);

        // add dM/dx \cdot \dot x
        MathLib::addToMatrix(Jac, {du - dv, 0.0, 0.0, dv - du});

        // add dK/dx \cdot \dot x
        MathLib::addToMatrix(Jac, {-du - dv, -du, 0.0, du + 2.0 * dv});

        // add -db/dx
        MathLib::addToMatrix(Jac,
                             {2.0 - 2.0 * v, -2.0 * u, 0.0, 2.0 * v - 5.0});

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

    bool requiresNormalization() const override { return false; }

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

const double ODETraits<ODE3>::t_end = std::numbers::pi / 2.;
// ODE 3 end //////////////////////////////////////////////////////
