// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseFlowWithPP
{
struct TwoPhaseFlowWithPPProcessData
{
    //! Specific body forces applied to solid and fluid.
    //! It is usually used to apply gravitational forces.
    //! A vector of displacement dimension's length.
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;

    //! Enables lumping of the mass matrix.
    bool const has_mass_lumping;
    ParameterLib::Parameter<double> const& temperature;
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
