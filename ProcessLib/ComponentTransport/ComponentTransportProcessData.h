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

#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
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
        MaterialLib::PorousMedium::PorousMediaProperties&&
            porous_media_properties_,
        ProcessLib::Parameter<double> const& fluid_reference_density_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperties>&&
            fluid_properties_,
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        ProcessLib::Parameter<double> const& retardation_factor_,
        ProcessLib::Parameter<double> const& decay_rate_,
        Eigen::VectorXd const& specific_body_force_,
        bool const has_gravity_)
        : porous_media_properties(std::move(porous_media_properties_)),
          fluid_reference_density(fluid_reference_density_),
          fluid_properties(std::move(fluid_properties_)),
          media_map(std::move(media_map_)),
          retardation_factor(retardation_factor_),
          decay_rate(decay_rate_),
          specific_body_force(specific_body_force_),
          has_gravity(has_gravity_)
    {
    }

    ComponentTransportProcessData(ComponentTransportProcessData&&) = default;

    //! Copies are forbidden.
    ComponentTransportProcessData(ComponentTransportProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(ComponentTransportProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ComponentTransportProcessData&&) = delete;

    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties;
    Parameter<double> const& fluid_reference_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> fluid_properties;
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    Parameter<double> const& retardation_factor;
    Parameter<double> const& decay_rate;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
