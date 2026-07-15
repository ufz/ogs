// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <string_view>

#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct Pipe;

/// Sample a scalar BHE parameter and reject non-finite or non-positive
/// values.  Plain `v <= 0` does NOT reject NaN (NaN comparisons are always
/// false) nor Inf (Inf <= 0 is false), so both slip through into std::log /
/// std::acosh and produce silent NaN in assembly.  Use this helper at every
/// BHE parameter evaluation site.
double sampleStrictPositive(ParameterLib::Parameter<double> const& param,
                            double const t,
                            ParameterLib::SpatialPosition const& pos,
                            std::string_view const param_role);

/// Assemble the shared "BHE geometry invalid" diagnostic prefix used by the
/// acosh-argument checks in the U-type resistance formulas. `equation` names
/// the formula, `cause` the physical degeneracy driving the argument to <= 1.
/// Kept in one place so the D/d0/distance message stays consistent across BHE
/// types.
std::string uTypeGeometryContext(ParameterLib::SpatialPosition const& pos,
                                 double const D,
                                 double const d0,
                                 double const distance,
                                 std::string_view const equation,
                                 std::string_view const cause);

/// Validate that `arg` is a valid input to std::acosh, i.e. finite and > 1.
/// The bound is strict: acosh(1) = 0 would make a thermal resistance zero and
/// produce an infinite (1/R) assembly coefficient.
/// On failure raises OGS_FATAL prefixed with the caller-assembled `context`,
/// which is expected to carry the element id and relevant geometry.
void checkAcoshArg(double const arg, std::string_view const context);

/// Validate that the borehole diameter D strictly exceeds `min_diameter`, the
/// lower bound the resistance formula requires (e.g. sqrt(2)*d0 for 1U, 2*d0
/// for 2U) so that std::log(D/min_diameter) stays real and positive.
void checkBoreholeVsPipeDiameter(double const D,
                                 double const min_diameter,
                                 ParameterLib::SpatialPosition const& pos,
                                 std::string_view const context);

/// Validate that a grout cross-sectional area (borehole area fraction minus
/// pipe outside area) is positive.  Returns the area on success, calls
/// OGS_FATAL on failure.
double checkedGroutArea(double const borehole_area_fraction,
                        double const pipe_outside_area,
                        ParameterLib::SpatialPosition const& pos);

/// Validate that the inlet and outlet pipes of a U-type BHE share the same
/// outside diameter. The Diersch (2011) grout and inter-grout resistance
/// formulas (R_g, R_ar) are written for a single pipe diameter d0, so a U-type
/// BHE must use one common outside diameter; the inlet and outlet may still
/// differ in wall thermal conductivity. Calls OGS_FATAL on mismatch. `context`
/// names the BHE type in the diagnostic message.
void checkEqualPipeOutsideDiameters(Pipe const& inlet, Pipe const& outlet,
                                    std::string_view const context);

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
