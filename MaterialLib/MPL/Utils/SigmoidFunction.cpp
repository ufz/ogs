/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 16, 2019, 3:40 PM
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
    return 1. / (1. + std::exp(k_ * (T - T_c_)));
}

double SigmoidFunction::dValue(double const& T) const
{
    double f = value(T);

    return -k_ * std::exp(k_ * (T - T_c_)) * (f * f);
}

double SigmoidFunction::d2Value(double const& T) const
{
    double f = value(T);
    double fT = dValue(T);

    return fT * (k_ + 2 * fT / f);
}
}  // namespace MaterialPropertyLib
