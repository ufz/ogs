/**
 * \file
 *
 * 2014/06/04 HS inital implementation
 * borehole heat exchanger abstract class
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
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
};
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
