/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TimeDiscretization.cpp
 *  Created on March 31, 2017, 10:30 AM
 */

#include "TimeDiscretization.h"
namespace NumLib
{
/**
 * Compute and return the relative change of solutions between two successive
 * time steps by \f$ e_n = \|u^{n+1}-u^{n}\|/\|u^{n+1}\| \f$.
 *
 * @param x      The current solution
 * @param x_old  The previous solution
 * @return       \f$ e_n = \|u^{n+1}-u^{n}\|/\|u^{n+1}\| \f$.
 *
 *  \warning the value of x_old is changed to x - x_old after this computation.
 */
static double computeRelativeError(GlobalVector const& x, GlobalVector& x_old,
                                   MathLib::VecNormType norm_type)
{
    if (norm_type == MathLib::VecNormType::INVALID)
        return 0.;

    const double norm2_x = MathLib::LinAlg::norm(x, norm_type);
    assert(norm2_x > std::numeric_limits<double>::epsilon());

    // dx = x - x_old --> x_old
    MathLib::LinAlg::axpy(x_old, -1.0, x);
    return MathLib::LinAlg::norm(x_old, norm_type) / norm2_x;
}

double BackwardEuler::getRelativeError(GlobalVector const& x,
                                       MathLib::VecNormType norm_type)
{
    return computeRelativeError(x, _x_old, norm_type);
}

double ForwardEuler::getRelativeError(GlobalVector const& x,
                                      MathLib::VecNormType norm_type)
{
    return computeRelativeError(x, _x_old, norm_type);
}

double CrankNicolson::getRelativeError(GlobalVector const& x,
                                       MathLib::VecNormType norm_type)
{
    return computeRelativeError(x, _x_old, norm_type);
}

double BackwardDifferentiationFormula::getRelativeError(
    GlobalVector const& x, MathLib::VecNormType norm_type)
{
    return (_xs_old.size() < _num_steps)
               ? 0.
               : computeRelativeError(x, *_xs_old[_offset], norm_type);
}

void BackwardDifferentiationFormula::pushState(const double,
                                               GlobalVector const& x,
                                               InternalMatrixStorage const&)
{
    namespace LinAlg = MathLib::LinAlg;
    // TODO use boost circular buffer?

    // until _xs_old is filled, lower-order BDF formulas are used.
    if (_xs_old.size() < _num_steps)
    {
        _xs_old.push_back(&NumLib::GlobalVectorProvider::provider.getVector(x));
    }
    else
    {
        LinAlg::copy(x, *_xs_old[_offset]);
        _offset =
            (_offset + 1) % _num_steps;  // treat _xs_old as a circular buffer
    }
}

namespace detail
{
//! Coefficients used in the backward differentiation formulas.
const double BDF_Coeffs[6][7] = {
    // leftmost column: weight of the solution at the new timestep
    // signs of columns > 1 are flipped compared to standard BDF tableaus
    {1.0, 1.0},
    {1.5, 2.0, -0.5},
    {11.0 / 6.0, 3.0, -1.5, 1.0 / 3.0},
    {25.0 / 12.0, 4.0, -3.0, 4.0 / 3.0, -0.25},
    {137.0 / 60.0, 5.0, -5.0, 10.0 / 3.0, -1.25, 1.0 / 5.0},
    {147.0 / 60.0, 6.0, -7.5, 20.0 / 3.0, -3.75, 6.0 / 5.0, -1.0 / 6.0}
    // coefficient of (for BDF(6), the oldest state, x_n, is always rightmost)
    //        x_+6, x_+5, x_+4,       x_+3,  x_+2, x_+1,     x_n
};
}

double BackwardDifferentiationFormula::getNewXWeight() const
{
    auto const k = eff_num_steps();
    return detail::BDF_Coeffs[k - 1][0] / _delta_t;
}

void BackwardDifferentiationFormula::getWeightedOldX(GlobalVector& y) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const k = eff_num_steps();
    auto const* const BDFk = detail::BDF_Coeffs[k - 1];

    // compute linear combination \sum_{i=0}^{k-1} BDFk_{k-i} \cdot x_{n+i}
    LinAlg::copy(*_xs_old[_offset], y);  // _xs_old[offset] = x_n
    LinAlg::scale(y, BDFk[k]);

    for (unsigned i = 1; i < k; ++i)
    {
        auto const off = (_offset + i) % k;
        LinAlg::axpy(y, BDFk[k - i], *_xs_old[off]);
    }

    LinAlg::scale(y, 1.0 / _delta_t);
}

}  // end of namespace NumLib