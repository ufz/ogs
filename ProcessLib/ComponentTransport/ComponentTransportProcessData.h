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

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/PorousMediaProperties.h"

namespace MaterialPropertyLib
{
class Medium;
}

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace ComponentTransport
{
struct ComponentTransportProcessData
{
    ComponentTransportProcessData(
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        ParameterLib::Parameter<double> const& decay_rate_,
        Eigen::VectorXd const& specific_body_force_, bool const has_gravity_,
        bool const non_advective_form_)
        : media_map(std::move(media_map_)),
          decay_rate(decay_rate_),
          specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          non_advective_form(non_advective_form_)
    {
    }

    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    ParameterLib::Parameter<double> const& decay_rate;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    bool const non_advective_form;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
