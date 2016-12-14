/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   ProcessType.h
 *
 * Created on December 14, 2016, 3:13 PM
 */

#ifndef OGS_PROCESS_TYPE_H
#define OGS_PROCESS_TYPE_H

namespace ProcessLib
{
enum class ProcessType
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
#endif /* OGS_PROCESS_TYPE_H */
