/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PorousMediaProperties.h"

#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
int PorousMediaProperties::getMaterialID(
    ParameterLib::SpatialPosition const& pos) const
{
    int const element_id = pos.getElementID().value();
    return _material_ids[element_id];
}

MaterialLib::PorousMedium::CapillaryPressureSaturation const&
PorousMediaProperties::getCapillaryPressureSaturationModel(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *_capillary_pressure_saturation_models[getMaterialID(pos)];
}
}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
