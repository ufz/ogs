/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MatrixTranslator.h"

#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& K,
             GlobalMatrix& A) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const dxdot_dx = _time_disc.getNewXWeight();

    // A = M * dxdot_dx + K
    LinAlg::copy(M, A);
    LinAlg::aypx(A, dxdot_dx, K);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto& tmp = NumLib::GlobalVectorProvider::provider.getVector(_tmp_id);
    _time_disc.getWeightedOldX(tmp);

    // rhs = M * weighted_old_x + b
    LinAlg::matMultAdd(M, tmp, b, rhs);

    NumLib::GlobalVectorProvider::provider.releaseVector(tmp);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

    // res = M * x_dot + K * x_curr - b
    LinAlg::matMult(M, xdot, res);  // the local vector x_dot seems to be
                                  // necessary because of this multiplication
    LinAlg::matMultAdd(K, x_curr, res, res);
    LinAlg::axpy(res, -1.0, b);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace LinAlg = MathLib::LinAlg;

    LinAlg::copy(Jac_in, Jac_out);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& /*K*/,
             GlobalMatrix& A) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const dxdot_dx = _fwd_euler.getNewXWeight();

    // A = M * dxdot_dx
    LinAlg::copy(M, A);
    LinAlg::scale(A, dxdot_dx);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& K,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto& tmp = NumLib::GlobalVectorProvider::provider.getVector(_tmp_id);
    _fwd_euler.getWeightedOldX(tmp);

    auto const& x_old = _fwd_euler.getXOld();

    // rhs = b + M * weighted_old_x - K * x_old
    LinAlg::matMult(K, x_old, rhs);        // rhs = K * x_old
    LinAlg::aypx(rhs, -1.0, b);            // rhs = b - K * x_old
    LinAlg::matMultAdd(M, tmp, rhs, rhs);  // rhs += M * weighted_old_x

    NumLib::GlobalVectorProvider::provider.releaseVector(tmp);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const& x_curr = _fwd_euler.getCurrentX(x_new_timestep);

    // res = M * x_dot + K * x_curr - b
    LinAlg::matMult(M, xdot, res);
    LinAlg::matMultAdd(K, x_curr, res, res);
    LinAlg::axpy(res, -1.0, b);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace LinAlg = MathLib::LinAlg;

    LinAlg::copy(Jac_in, Jac_out);
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& K,
             GlobalMatrix& A) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const dxdot_dx = _crank_nicolson.getNewXWeight();
    auto const theta = _crank_nicolson.getTheta();

    // A = theta * (M * dxdot_dx + K) + dxdot_dx * _M_bar
    LinAlg::copy(M, A);
    LinAlg::aypx(A, dxdot_dx, K);  // A = M * dxdot_dx + K

    LinAlg::scale(A, theta);            // A *= theta
    LinAlg::axpy(A, dxdot_dx, _M_bar);  // A += dxdot_dx * _M_bar
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto& tmp = NumLib::GlobalVectorProvider::provider.getVector(_tmp_id);
    _crank_nicolson.getWeightedOldX(tmp);

    auto const theta = _crank_nicolson.getTheta();

    // rhs = theta * (b + M * weighted_old_x) + _M_bar * weighted_old_x -
    // _b_bar;
    LinAlg::matMultAdd(M, tmp, b, rhs);  // rhs = b + M * weighted_old_x

    LinAlg::scale(rhs, theta);                  // rhs *= theta
    LinAlg::matMultAdd(_M_bar, tmp, rhs, rhs);  // rhs += _M_bar * weighted_old_x
    LinAlg::axpy(rhs, -1.0, _b_bar);            // rhs -= b

    NumLib::GlobalVectorProvider::provider.releaseVector(tmp);
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const& x_curr = _crank_nicolson.getCurrentX(x_new_timestep);
    auto const theta = _crank_nicolson.getTheta();

    // res = theta * (M * x_dot + K*x_curr - b) + _M_bar * x_dot + _b_bar
    LinAlg::matMult(M, xdot, res);            // res = M * x_dot
    LinAlg::matMultAdd(K, x_curr, res, res);  // res += K * x_curr
    LinAlg::axpy(res, -1.0, b);               // res = M * x_dot + K * x_curr - b

    LinAlg::aypx(res, theta, _b_bar);            // res = res * theta + _b_bar
    LinAlg::matMultAdd(_M_bar, xdot, res, res);  // rs += _M_bar * x_dot
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const dxdot_dx = _crank_nicolson.getNewXWeight();
    auto const theta = _crank_nicolson.getTheta();

    // J = theta * Jac + dxdot_dx * _M_bar
    LinAlg::copy(Jac_in, Jac_out);
    LinAlg::scale(Jac_out, theta);
    LinAlg::axpy(Jac_out, dxdot_dx, _M_bar);
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    pushMatrices(GlobalMatrix const& M, GlobalMatrix const& K,
                 GlobalVector const& b)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const theta = _crank_nicolson.getTheta();

    // Note: using x_old here is correct, since this method is called from
    // within
    //       CrankNicolson::pushState() __after__ x_old has been updated to the
    //       result
    //       from the timestep just finished.
    auto const& x_old = _crank_nicolson.getXOld();

    // _M_bar = (1.0-theta) * M;
    LinAlg::copy(M, _M_bar);
    LinAlg::scale(_M_bar, 1.0 - theta);

    // _b_bar = (1.0-theta) * (K * x_old - b)
    LinAlg::matMult(K, x_old, _b_bar);
    LinAlg::axpy(_b_bar, -1.0, b);
    LinAlg::scale(_b_bar, 1.0 - theta);
}

}  // NumLib
