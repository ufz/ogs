// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
/// Sample a scalar BHE parameter and reject non-finite or non-positive
/// values.  Plain `v <= 0` does NOT reject NaN (NaN comparisons are always
/// false) nor Inf (Inf <= 0 is false), so both slip through into std::log /
/// std::acosh and produce silent NaN in assembly.  Use this helper at every
/// BHE parameter evaluation site.
double sampleStrictPositive(ParameterLib::Parameter<double> const& param,
                            double const t,
                            ParameterLib::SpatialPosition const& pos,
                            std::string_view const param_role);

/// Validate that `arg` is a valid input to std::acosh, i.e. finite and >= 1.
/// On failure raises OGS_FATAL prefixed with the caller-assembled `context`,
/// which is expected to carry the element id and relevant geometry.
void checkAcoshArg(double const arg, std::string_view const context);

/// Validate that the borehole diameter D strictly exceeds d0 (so std::log(D/d0)
/// stays real and positive).
void checkBoreholeVsPipeDiameter(double const D,
                                 double const d0,
                                 ParameterLib::SpatialPosition const& pos,
                                 std::string_view const context);

/// Validate that a grout cross-sectional area (borehole area fraction minus
/// pipe outside area) is positive.  Returns the area on success, calls
/// OGS_FATAL on failure.
double checkedGroutArea(double const borehole_area_fraction,
                        double const pipe_outside_area,
                        ParameterLib::SpatialPosition const& pos);

/// Validate a scalar, time-invariant BHE parameter (borehole diameter, pipe
/// wall thermal conductivity, ...) at config time. Rejects genuinely
/// time-varying parameter types (and warns for FunctionParameter, which is only
/// usable when spatial-only), rejects MeshNodeParameter, and requires a single
/// global component. Calls OGS_FATAL on violation. `role` names the parameter
/// in the diagnostic messages.
void validateScalarTimeInvariantBHEParameter(
    ParameterLib::Parameter<double> const& param, std::string_view const role);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
