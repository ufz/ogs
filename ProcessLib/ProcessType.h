/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   ProcessType.h
 *
 * Created on December 14, 2016, 3:13 PM
 */


#pragma once

namespace ProcessLib
{
enum ProcessType
{
    GroundwaterFlowProcess,
    LiquidFlowProcess,
    RichardsFlowProcess,
    TwoPhaseFlowWithPPProcess,
    HeatConductionProcess,
    HTProcess,
    HydroMechanicsProcess,
    SmallDeformationProcess,
    TESProcess
};
}
