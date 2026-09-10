// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MaterialPropertyLib
{
/// Steam dryness \f$x\f$, the vapour mass fraction, from the specific enthalpy
/// of the mixture and the saturation enthalpies of the two phases.
///
/// The dryness is capped at the two ends of the two-phase range, but the two
/// caps do not carry the same weight.
///
/// Below the range, \f$h < h_l\f$, the state is subcooled liquid. The callers
/// detect it as a dryness of zero and then re-evaluate the liquid density and
/// the temperature off the saturation line, so the cap only removes a negative
/// vapour mass fraction from an otherwise fully represented state.
///
/// Above the range, \f$h > h_v\f$, the state is superheated steam and there is
/// no counterpart: the mixture is described by the saturation line alone, so
/// such a section keeps the saturation temperature and the saturation vapour
/// density, and its superheat is lost rather than represented. The cap is
/// still applied, because the drift-flux closure is only defined on \f$[0,
/// 1]\f$, but it is reported so that it does not pass silently.
///
/// A non-finite ratio, which the saturation enthalpies produce where they
/// coincide at the critical point, passes through both caps unchanged and is
/// rejected by the closure that consumes the dryness.
double steamDryness(double const enthalpy, double const h_sat_liquid,
                    double const h_sat_vapour);
}  // namespace MaterialPropertyLib
