// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BHEParameterValidation.h"

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <cmath>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "ParameterLib/CurveScaledParameter.h"
#include "ParameterLib/MeshNodeParameter.h"
#include "ParameterLib/TimeDependentHeterogeneousParameter.h"
#include "Pipe.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
double sampleStrictPositive(ParameterLib::Parameter<double> const& param,
                            double const t,
                            ParameterLib::SpatialPosition const& pos,
                            std::string_view const param_role)
{
    double const v = param(t, pos)[0];
    if (!std::isfinite(v) || v <= 0.0)
    {
        OGS_FATAL(
            "BHE {:s} ('{}') is non-finite or non-positive (value={:g}) at "
            "{}.",
            param_role, param.name, v, pos);
    }
    return v;
}

std::string uTypeGeometryContext(ParameterLib::SpatialPosition const& pos,
                                 double const D,
                                 double const d0,
                                 double const distance,
                                 std::string_view const equation,
                                 std::string_view const cause)
{
    return fmt::format(
        "BHE geometry invalid at {} ({:s}): D={:g}, d0={:g}, distance={:g}. "
        "{:s}",
        pos, equation, D, d0, distance, cause);
}

void checkAcoshArg(double const arg, std::string_view const context)
{
    if (!std::isfinite(arg) || arg <= 1.0)
    {
        OGS_FATAL("{:s} -> acosh argument = {:g} (must be finite and > 1).",
                  context, arg);
    }
}

void checkBoreholeVsPipeDiameter(double const D,
                                 double const min_diameter,
                                 ParameterLib::SpatialPosition const& pos,
                                 std::string_view const context)
{
    if (!(D > min_diameter) || !std::isfinite(D / min_diameter))
    {
        OGS_FATAL(
            "BHE geometry invalid at {} ({:s}): borehole diameter "
            "D={:g} must strictly exceed the minimum diameter {:g} required by "
            "this formula.",
            pos, context, D, min_diameter);
    }
}

double checkedGroutArea(double const borehole_area_fraction,
                        double const pipe_outside_area,
                        ParameterLib::SpatialPosition const& pos)
{
    double const grout_area = borehole_area_fraction - pipe_outside_area;
    if (grout_area <= 0)
    {
        OGS_FATAL(
            "Non-positive grout cross-sectional area at {}. "
            "Borehole diameter is too small for the pipe dimensions.",
            pos);
    }
    return grout_area;
}

void checkEqualPipeOutsideDiameters(Pipe const& inlet, Pipe const& outlet,
                                    std::string_view const context)
{
    double const d_inlet = inlet.outsideDiameter();
    double const d_outlet = outlet.outsideDiameter();
    // The diameters are user-supplied constants, so allow for round-off but
    // reject genuinely different pipes.
    if (std::abs(d_inlet - d_outlet) >
        1e-12 * std::max(std::abs(d_inlet), std::abs(d_outlet)))
    {
        OGS_FATAL(
            "{:s}: inlet and outlet pipe outside diameters differ "
            "(inlet={:g}, outlet={:g}). U-type BHE thermal resistance formulas "
            "assume a single pipe diameter; inlet and outlet may differ only "
            "in wall thermal conductivity.",
            context, d_inlet, d_outlet);
    }
}

void validateScalarTimeInvariantBHEParameter(
    ParameterLib::Parameter<double> const& param, std::string_view const role)
{
    // Borehole geometry properties are sampled once at t=0 and are physically
    // time-invariant, so reject parameter types that genuinely vary in time.
    // Parameter::isTimeDependent() is too conservative to use as the sole
    // criterion here: a FunctionParameter reports true even for a purely
    // spatial expression. We therefore deny-list only the parameter types that
    // are always genuinely time-varying, and merely warn (below) for anything
    // else that reports time-dependence.
    if (dynamic_cast<ParameterLib::CurveScaledParameter<double> const*>(
            &param) != nullptr ||
        dynamic_cast<ParameterLib::TimeDependentHeterogeneousParameter const*>(
            &param) != nullptr)
    {
        OGS_FATAL(
            "BHE {:s} '{}' uses a time-varying parameter type, but borehole "
            "geometry properties are sampled once at t=0 and are physically "
            "time-invariant. Use ConstantParameter, MeshElementParameter, or "
            "a spatial-only FunctionParameter.",
            role, param.name);
    }
    if (param.isTimeDependent())
    {
        WARN(
            "BHE {:s} '{}' reports time-dependence. This is expected for a "
            "FunctionParameter, whose time-dependence flag is set even for a "
            "purely spatial expression, so a spatial-only expression is fine. "
            "Only the value at t=0 is used; ensure the expression does not "
            "actually depend on time.",
            role, param.name);
    }

    if (dynamic_cast<ParameterLib::MeshNodeParameter<double> const*>(&param))
    {
        OGS_FATAL(
            "BHE parameter '{}' must not be a MeshNodeParameter. "
            "Use MeshElementParameter for spatially varying BHE properties.",
            param.name);
    }

    if (param.getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "BHE {:s} '{}' must be scalar (1 component), got {:d} components.",
            role, param.name, param.getNumberOfGlobalComponents());
    }
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
