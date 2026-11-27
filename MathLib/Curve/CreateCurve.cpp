// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateCurve.h"

#include <string>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "Curve.h"
#include "PiecewiseLinearCurveConfig.h"

namespace MathLib
{
std::unique_ptr<Curve> createCurve(std::string_view const name,
                                   BaseLib::ConfigTree const& config)
{
    auto const interpolation =
        //! \ogs_file_param{prj__curves__curve__interpolation}
        config.getConfigParameter<std::string>("interpolation", "linear");

    auto const left_end_slope =
        //! \ogs_file_param{prj__curves__curve__left_end_slope}
        config.getConfigParameterOptional<double>("left_end_slope");
    auto const right_end_slope =
        //! \ogs_file_param{prj__curves__curve__right_end_slope}
        config.getConfigParameterOptional<double>("right_end_slope");

    if (interpolation == "linear" && (left_end_slope || right_end_slope))
    {
        OGS_FATAL(
            "The end slopes of curve '{:s}' can only be given for the "
            "interpolation 'pchip'.",
            name);
    }
    if (interpolation != "linear" && interpolation != "pchip")
    {
        OGS_FATAL(
            "Unknown interpolation '{:s}' for curve '{:s}'. Available are "
            "'linear' and 'pchip'.",
            interpolation, name);
    }

    auto [xs, ys] = parsePiecewiseLinearCurveConfig(config);
    if (interpolation == "linear")
    {
        return std::make_unique<Curve>(PiecewiseLinear{}, name, std::move(xs),
                                       std::move(ys));
    }
    return std::make_unique<Curve>(Pchip{.left_end_slope = left_end_slope,
                                         .right_end_slope = right_end_slope},
                                   name, std::move(xs), std::move(ys));
}
}  // namespace MathLib
