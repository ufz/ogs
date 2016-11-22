/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreatePiecewiseLinearCurve.h
 *
 * Created on November 11, 2016, 10:49 AM
 */

#ifndef OGS_CREATE_PIECEWISE_LINEAR_CURVE_H
#define OGS_CREATE_PIECEWISE_LINEAR_CURVE_H

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
class PiecewiseLinearCurve;

/** Create a curve
 *  \param config             ConfigTree object has a tag of <curve>
 *  \param check_monotonicity A flag for checking the monotonicity of the curve.
 */
std::unique_ptr<PiecewiseLinearCurve> createPiecewiseLinearCurve(
    BaseLib::ConfigTree const& config, bool const check_monotonicity = false);
};

#endif /* OGS_CREATE_PIECEWISE_LINEAR_CURVE_H */
