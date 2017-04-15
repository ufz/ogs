/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#pragma once

namespace ProcessLib
{

/**
 * y ... variable in global matrix
 * x ... "physical" process variable in local assembly
 *
 * x = exp(y), dx = dx/dy * dy
 * dx/dy = exp(y) = x
 */
struct TrafoLog
{
    static const bool constrained = true;

    /// Converts global matrix entry to "physical" variable
    /// used in local assembly.
    static double x(const double y) { return std::exp(y); }

    /// Derivative of the "physical" variable x w.r.t. y.
    /// the argument is x!
    static double dxdy(const double x) { return x; }
};

struct TrafoIdentity
{
    static const bool constrained = false;

    /// Converts global matrix entry to "physical" variable
    /// used in local assembly.
    static double x(const double y) { return y; }

    /// Derivative of the "physical" variable x w.r.t. y.
    /// the argument is x!
    static double dxdy(const double /*x*/) { return 1.0; }
};

struct TrafoTanh
{
    static const bool constrained = true;

    /// Converts global matrix entry to "physical" variable
    /// used in local assembly.
    static double x(const double y) { return 0.5 * std::tanh(y) + 0.5; }

    /// Derivative of the "physical" variable x w.r.t. y.
    /// the argument is x!
    static double dxdy(const double x) { return 2.0*x*(1.0-x); }
};

struct TrafoScale
{
    TrafoScale(const double factor) : _factor{factor} {}

    static const bool constrained = false;

    /// Converts global matrix entry to "physical" variable
    /// used in local assembly.
    double x(const double y) const { return y * _factor; }

    double y(const double x) const { return x / _factor; }

    /// Derivative of the "physical" variable x w.r.t. y.
    /// the argument is x!
    double dxdy(const double /*x*/) const { return _factor; }

private:
    double _factor;
};

using Trafo = ProcessLib::TrafoScale;
}
