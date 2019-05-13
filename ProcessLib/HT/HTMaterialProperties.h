/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace HT
{
struct HTMaterialProperties final
{
    HTMaterialProperties(
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        bool const has_fluid_thermal_expansion_,
        ParameterLib::Parameter<double> const& solid_thermal_expansion_,
        ParameterLib::Parameter<double> const& biot_constant_,
        Eigen::VectorXd specific_body_force_,
        bool const has_gravity_)
        : media_map(std::move(media_map_)),
          has_fluid_thermal_expansion(has_fluid_thermal_expansion_),
          solid_thermal_expansion(solid_thermal_expansion_),
          biot_constant(biot_constant_),
          specific_body_force(std::move(specific_body_force_)),
          has_gravity(has_gravity_)
    {
    }

    HTMaterialProperties(HTMaterialProperties&&) = delete;
    HTMaterialProperties(HTMaterialProperties const&) = delete;
    void operator=(HTMaterialProperties&&) = delete;
    void operator=(HTMaterialProperties const&) = delete;

    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;

    bool const has_fluid_thermal_expansion;
    ParameterLib::Parameter<double> const& solid_thermal_expansion;
    ParameterLib::Parameter<double> const& biot_constant;

    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace HT
}  // namespace ProcessLib
