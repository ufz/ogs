/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MatrixTranslator.h"

#include "MathLib/LinAlg/BLAS.h"

namespace NumLib
{
void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& K,
             GlobalMatrix& A) const
{
    namespace BLAS = MathLib::BLAS;

    auto const dxdot_dx = _time_disc.getNewXWeight();

    // A = M * dxdot_dx + K
    BLAS::copy(M, A);
    BLAS::aypx(A, dxdot_dx, K);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace BLAS = MathLib::BLAS;

    auto& tmp = MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(_tmp_id);
    _time_disc.getWeightedOldX(tmp);

    // rhs = M * weighted_old_x + b
    BLAS::matMultAdd(M, tmp, b, rhs);

    MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(tmp);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace BLAS = MathLib::BLAS;

    auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

    // res = M * x_dot + K * x_curr - b
    BLAS::matMult(M, xdot, res);  // the local vector x_dot seems to be
                                  // necessary because of this multiplication
    BLAS::matMultAdd(K, x_curr, res, res);
    BLAS::axpy(res, -1.0, b);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace BLAS = MathLib::BLAS;

    BLAS::copy(Jac_in, Jac_out);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& /*K*/,
             GlobalMatrix& A) const
{
    namespace BLAS = MathLib::BLAS;

    auto const dxdot_dx = _fwd_euler.getNewXWeight();

    // A = M * dxdot_dx
    BLAS::copy(M, A);
    BLAS::scale(A, dxdot_dx);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& K,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace BLAS = MathLib::BLAS;

    auto& tmp = MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(_tmp_id);
    _fwd_euler.getWeightedOldX(tmp);

    auto const& x_old = _fwd_euler.getXOld();

    // rhs = b + M * weighted_old_x - K * x_old
    BLAS::matMult(K, x_old, rhs);        // rhs = K * x_old
    BLAS::aypx(rhs, -1.0, b);            // rhs = b - K * x_old
    BLAS::matMultAdd(M, tmp, rhs, rhs);  // rhs += M * weighted_old_x

    MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(tmp);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace BLAS = MathLib::BLAS;

    auto const& x_curr = _fwd_euler.getCurrentX(x_new_timestep);

    // res = M * x_dot + K * x_curr - b
    BLAS::matMult(M, xdot, res);
    BLAS::matMultAdd(K, x_curr, res, res);
    BLAS::axpy(res, -1.0, b);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace BLAS = MathLib::BLAS;

    BLAS::copy(Jac_in, Jac_out);
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& K,
             GlobalMatrix& A) const
{
    namespace BLAS = MathLib::BLAS;

    auto const dxdot_dx = _crank_nicolson.getNewXWeight();
    auto const theta = _crank_nicolson.getTheta();

    // A = theta * (M * dxdot_dx + K) + dxdot_dx * _M_bar
    BLAS::copy(M, A);
    BLAS::aypx(A, dxdot_dx, K);  // A = M * dxdot_dx + K

    BLAS::scale(A, theta);            // A *= theta
    BLAS::axpy(A, dxdot_dx, _M_bar);  // A += dxdot_dx * _M_bar
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace BLAS = MathLib::BLAS;

    auto& tmp = MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(_tmp_id);
    _crank_nicolson.getWeightedOldX(tmp);

    auto const theta = _crank_nicolson.getTheta();

    // rhs = theta * (b + M * weighted_old_x) + _M_bar * weighted_old_x -
    // _b_bar;
    BLAS::matMultAdd(M, tmp, b, rhs);  // rhs = b + M * weighted_old_x

    BLAS::scale(rhs, theta);                  // rhs *= theta
    BLAS::matMultAdd(_M_bar, tmp, rhs, rhs);  // rhs += _M_bar * weighted_old_x
    BLAS::axpy(rhs, -1.0, _b_bar);            // rhs -= b

    MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(tmp);
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace BLAS = MathLib::BLAS;

    auto const& x_curr = _crank_nicolson.getCurrentX(x_new_timestep);
    auto const theta = _crank_nicolson.getTheta();

    // res = theta * (M * x_dot + K*x_curr - b) + _M_bar * x_dot + _b_bar
    BLAS::matMult(M, xdot, res);            // res = M * x_dot
    BLAS::matMultAdd(K, x_curr, res, res);  // res += K * x_curr
    BLAS::axpy(res, -1.0, b);               // res = M * x_dot + K * x_curr - b

    BLAS::aypx(res, theta, _b_bar);            // res = res * theta + _b_bar
    BLAS::matMultAdd(_M_bar, xdot, res, res);  // rs += _M_bar * x_dot
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace BLAS = MathLib::BLAS;

    auto const dxdot_dx = _crank_nicolson.getNewXWeight();
    auto const theta = _crank_nicolson.getTheta();

    // J = theta * Jac + dxdot_dx * _M_bar
    BLAS::copy(Jac_in, Jac_out);
    BLAS::scale(Jac_out, theta);
    BLAS::axpy(Jac_out, dxdot_dx, _M_bar);
}

void MatrixTranslatorCrankNicolson<
    ODESystemTag::FirstOrderImplicitQuasilinear>::
    pushMatrices(GlobalMatrix const& M, GlobalMatrix const& K,
                 GlobalVector const& b)
{
    namespace BLAS = MathLib::BLAS;

    auto const theta = _crank_nicolson.getTheta();

    // Note: using x_old here is correct, since this method is called from
    // within
    //       CrankNicolson::pushState() __after__ x_old has been updated to the
    //       result
    //       from the timestep just finished.
    auto const& x_old = _crank_nicolson.getXOld();

    // _M_bar = (1.0-theta) * M;
    BLAS::copy(M, _M_bar);
    BLAS::scale(_M_bar, 1.0 - theta);

    // _b_bar = (1.0-theta) * (K * x_old - b)
    BLAS::matMult(K, x_old, _b_bar);
    BLAS::axpy(_b_bar, -1.0, b);
    BLAS::scale(_b_bar, 1.0 - theta);
}

}  // NumLib
