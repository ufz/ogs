/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "PorousMediaProperties.h"

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;
}

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace RichardsComponentTransport
{
struct RichardsComponentTransportProcessData
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    PorousMediaProperties porous_media_properties;
    ParameterLib::Parameter<double> const& fluid_reference_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> fluid_properties;
    ParameterLib::Parameter<double> const& molecular_diffusion_coefficient;
    ParameterLib::Parameter<double> const& solute_dispersivity_longitudinal;
    ParameterLib::Parameter<double> const& solute_dispersivity_transverse;
    ParameterLib::Parameter<double> const& retardation_factor;
    ParameterLib::Parameter<double> const& decay_rate;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
