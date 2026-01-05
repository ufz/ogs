// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 *
 * 1) Diersch_2011_CG
 * Two very important references to understand this class implementations are:
 * Diersch, H-JG, D. Bauer, W. Heidemann, Wolfram Rühaak, and Peter Schätzl.
 * Finite element modeling of borehole heat exchanger systems:
 * Part 1. Fundamentals, Computers & Geosciences,
 * Volume 37, Issue 8, August 2011, Pages 1122-1135, ISSN 0098-3004,
 * http://dx.doi.org/10.1016/j.cageo.2010.08.003.
 *
 * 2) FEFLOW_2014_Springer
 * FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport in Porous
 * and Fractured Media Diersch, Hans-Joerg, 2014, XXXV, 996 p, Springer.
 */

#pragma once

#include "BoreholeGeometry.h"
#include "FlowAndTemperatureControl.h"
#include "GroutParameters.h"
#include "RefrigerantProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct BHECommon
{
    BoreholeGeometry const borehole_geometry;
    RefrigerantProperties const refrigerant;
    GroutParameters const grout;
    FlowAndTemperatureControl const flowAndTemperatureControl;
    bool const use_python_bcs;
    constexpr bool isPowerBC() const
    {
        return std::visit([](auto const& ftc) { return ftc.is_power_bc; },
                          flowAndTemperatureControl);
    }
};
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
