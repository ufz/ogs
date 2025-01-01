/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SpecificBodyForceData
{
    Eigen::Vector<double, DisplacementDim> specific_body_force;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
