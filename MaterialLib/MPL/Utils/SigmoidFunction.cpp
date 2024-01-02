/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 20, 2022
 */

#include "SigmoidFunction.h"

namespace MaterialPropertyLib
{
SigmoidFunction::SigmoidFunction(double const k, double const T_c)
    : k_(k), T_c_(T_c)
{
}

double SigmoidFunction::value(double const& T) const
{
    double const x = k_ * (T - T_c_);

    // Cutting off at the last x producing a (non-normal) return value because
    // std::exp(x) produces +infinity beyond this point, approximately 709.78.
    // The reason for using hex representation is that the value is unambiguous.
    if (x > 0x1.62e42fefa39efp+9)
    {
        return 0.;
    }

    return 1. / (1. + std::exp(x));
}

double SigmoidFunction::dValue(double const& T) const
{
    double const f = value(T);

    // Skip any further computations because of f^2 below.
    if (f * f == 0)
    {
        return 0;
    }

    return -k_ * std::exp(k_ * (T - T_c_)) * (f * f);
}

double SigmoidFunction::d2Value(double const& T) const
{
    double const fT = dValue(T);

    if (fT == 0)
    {
        return 0;
    }

    double const f = value(T);
    return fT * (k_ + 2 * fT / f);
}
}  // namespace MaterialPropertyLib
