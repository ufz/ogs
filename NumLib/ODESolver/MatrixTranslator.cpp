/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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

    double const dt = _time_disc.getCurrentTimeIncrement();

    // A = M * 1/dt + K
    LinAlg::copy(K, A);
    LinAlg::axpy(A, 1. / dt, M);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
               const GlobalVector& b, const GlobalVector& x_prev,
               GlobalVector& rhs) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto& tmp = NumLib::GlobalVectorProvider::provider.getVector(_tmp_id);
    _time_disc.getWeightedOldX(tmp, x_prev);

    // rhs = M * weighted_old_x + b
    LinAlg::matMultAdd(M, tmp, b, rhs);

    NumLib::GlobalVectorProvider::provider.releaseVector(tmp);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_curr,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace LinAlg = MathLib::LinAlg;

    // res = M * x_dot + K * x_curr - b
    LinAlg::copy(b, res);
    LinAlg::scale(res, -1.0);

    LinAlg::matMultAdd(M, xdot, res,
                       res);  // the local vector x_dot seems to be
                              // necessary because of this multiplication
    LinAlg::matMultAdd(K, x_curr, res, res);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace LinAlg = MathLib::LinAlg;

    LinAlg::copy(Jac_in, Jac_out);
}

}  // namespace NumLib
