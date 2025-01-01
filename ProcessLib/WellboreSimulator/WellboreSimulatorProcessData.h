/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
#include <memory>
#include <utility>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ParameterLib/Parameter.h"
#include "ReservoirProperties.h"
#include "WellboreGeometry.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
struct WellboreSimulatorProcessData final
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
    Eigen::VectorXd const specific_body_force;
    WellboreGeometry wellbore;
    ParameterLib::Parameter<double> const& well_ref_pressure;
    ParameterLib::Parameter<double> const& well_ref_enthalpy;
    ReservoirProperties reservoir_properties;
    ParameterLib::Parameter<double> const& productivity_index;

    bool const has_heat_exchange_with_formation;
    bool const has_gravity;

    MeshLib::PropertyVector<double>* mesh_prop_density = nullptr;
};

}  // namespace WellboreSimulator
}  // namespace ProcessLib
