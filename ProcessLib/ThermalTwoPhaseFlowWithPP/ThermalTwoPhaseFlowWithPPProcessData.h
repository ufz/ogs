// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace ThermalTwoPhaseFlowWithPP
{
struct ThermalTwoPhaseFlowWithPPProcessData
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;
    bool const has_mass_lumping;
};

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
